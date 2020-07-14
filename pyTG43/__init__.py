from glob import glob
import numpy as np
import pydicom
from pydicom.tag import Tag
from terminaltables import SingleTable
from matplotlib.path import Path
import xlrd, math, os
from multiprocessing import Pool
import itertools


def tpsComp(rp, rs, directory):
    """Calculate and compare dose at reference points with TPS.

    Args:
        rp: pydicom object of plan (RP) file.
        rs: pydicom object of structure (RS) file.
        directory: directory containing source data.

    Returns:
        points: list of DosePoint objects.
    """

    points = []

    source = Source(rp, directory)
    plan = Plan(source, rp, rs)

    for p in rp[0x300A, 0x10]:
        if Tag(0x300A, 0x18) in p.keys():
            x, y, z = p[0x300A, 0x18].value
            name = p[0x300A, 0x16].value
            ref = p[0x300A, 0x12].value
            points.append(DosePoint([x / 10, y / 10, z / 10], source, plan, name, ref))

    table_data = [["Name", "X", "Y", "Z", "TPS (Gy)", "Calc (Gy)", "Diff (%)"]]

    for pt in points:
        table_data.append(
            [
                pt.name,
                "{:.2f}".format(pt.x),
                "{:.2f}".format(pt.y),
                "{:.2f}".format(pt.z),
                "{:.3f}".format(pt.tpsdose),
                "{:.3f}".format(pt.dose),
                "{:.3f}".format((1 - (pt.dose / pt.tpsdose)) * 100),
            ]
        )

    print(SingleTable(table_data).table)
    return points


def pcalc(pt):
    """ Parallel DosePoint calculation helper function. """
    return DosePoint(pt, source, plan).dose


def calcDVHs(sourcei, plani, maxd, names):
    """Calculate cumulative DVHs for structures.

    Atributes:
        sourcei: pyTG43.Source object.
        plani: pyTG43.Plan object.
        maxd: max value for DVH calculation.
        names: list of names of structures to calculate DVH for.

    Returns:
        roi.dvh: array of cumulative DVH for each ROI in names.
    """
    global source
    global plan

    source = sourcei
    plan = plani

    if os.name != "nt":
        pool = Pool()
    for roi in plan.ROIs:
        if roi.name.lower() in [x.lower() for x in names]:
            if os.name == "nt":
                dvh = [DosePoint(pt, source, plan).dose for pt in roi.dvhpts]
            else:
                dvh = pool.map(pcalc, roi.dvhpts)
            n, bins = np.histogram(dvh, 100, range=(0, maxd))
            dvh = np.cumsum(n[::-1])[::-1]
            dvh = dvh / dvh.max() * 100
            roi.dvh = np.column_stack((bins[:-1], dvh))

    if os.name != "nt":
        pool.close()
        pool.join()


class Source(object):
    """Source parameters object.

    Attributes:
        L: length (cm).
        Delta: dose rate constant (cGy h-1 U-1).
        Sk: RAKR.
        Fi(r,theta): anisotropy user data interpolation function.
        gi(r): radial dose user data interpolation function.
        G0: reference geometry function.
    """

    def __init__(self, rp, directory):
        """
        Args:
            rp: pydicom object of plan (RP) file.
            directory: directory containing source data.
        """
        r0 = 1
        theta0 = np.pi / 2

        self.Sk = rp.SourceSequence[0].ReferenceAirKermaRate

        fname = glob(directory + "/*" + rp.BrachyTreatmentType.lower() + "*.xls")[0]
        wb = xlrd.open_workbook(fname)
        sh = wb.sheets()[-1]

        self.L = sh.row(9)[2].value
        self.Delta = sh.row(4)[2].value

        col = 5
        while sh.row(10)[col].ctype == 2:
            col += 1

        Fi = np.ones((1, col - 5))
        gi = np.ones((1, 2))
        Ft = []
        Fr = [r.value for r in sh.row(10)[5:col]]

        for row in np.arange(11, sh.nrows):
            if sh.row(row)[1].ctype != 0:
                gi = np.vstack([gi, np.array([r.value for r in sh.row(row)[1:3]])])
            Fi = np.vstack([Fi, np.array([r.value for r in sh.row(row)[5:col]])])
            Ft.append(sh.row(row)[4].value)

        Fi = Fi[1:, :]
        gi = gi[1:, :]

        self.Fi = bilinearinterp(Fr, Ft, Fi)
        self.gi = fastinterp(gi[:, 0], gi[:, 1])

        self.G0 = 2 * np.arctan((self.L / 2) / r0) / (self.L * r0 * np.sin(theta0))

    def F(self, r, theta):
        """Anisotropy function.

        Args:
            r: distance (cm).
            theta: angle (radians).
        """
        if r > 10:
            return self.Fi(10, theta)
        else:
            return self.Fi(r, theta)

    def g(self, r):
        """Radial dose function.

        Args:
            r: distance (cm).
        """
        if r < 0.15:
            return self.gi(0.15)
        elif r > 10:
            return self.gi(8) * np.exp(
                (r - 8) / (10 - 8) * (np.log(self.gi(10)) - np.log(self.gi(8)))
            )
        else:
            return self.gi(r)


class Dwell(object):
    """Planned dwell point.

    Attributes:
        x, y, z: x, y and z coordinates (cm).
        middle: (x,y,z) coordinates of middle of source (cm).
        t: dwell time (s).
        ends: list containing coordinates of each endd of the source (cm).
        rotation: dwell position orientation angles for egs_brachy.
    """

    def __init__(self, coords, t, L, app):
        """
        Args:
            coords: list containing x, y, and z coordinates of dwell point.
            t: dwell time.
            L: source length.
            app: applicator ROI object dwell point is a part of.
        """

        self.x, self.y, self.z = coords
        self.middle = np.array(coords)
        self.t = t
        self.get_source_position(app, L)

    def get_source_position(self, app, L):
        """Calculates coordinates of the ends of the source.

        Args:
            app: applicator ROI object dwell point is a part of.
            L: source length (cm).
        """
        smallest = +np.inf
        for i, p in enumerate(app.oldcoords(app.coords) / 10):
            if i < len(app.coords) - 1:
                q = app.oldcoords(app.coords)[i + 1] / 10
                online = (
                    np.round(euclidzip(p, self.middle) + euclidzip(q, self.middle), 4)
                ) - np.round(euclidzip(p, q), 4)
                if online < smallest:
                    smallest = (
                        np.round(
                            euclidzip(p, self.middle) + euclidzip(q, self.middle), 4
                        )
                    ) - np.round(euclidzip(p, q), 4)
                    closest = [p, q]

        v = closest[0] - closest[1]
        r = euclidzip(closest[0], closest[1])
        nv = v / r

        cosines = [np.arccos(x) for x in nv]
        alpha = np.arccos(
            (np.cos(cosines[1])) / (np.cos(np.arcsin(-np.cos(cosines[0]))))
        )
        beta = np.arcsin(-np.cos(cosines[0]))

        if abs(np.cos(cosines[2]) - np.sin(alpha) * np.cos(beta)) > 1e-15:
            v = closest[1] - closest[0]
            r = euclidzip(closest[1], closest[0])
            nv = v / r
            cosines = [np.arccos(x) for x in nv]
            alpha = (
                np.arccos(
                    (np.cos(cosines[1])) / (np.cos(np.arcsin(-np.cos(cosines[0]))))
                )
                + np.pi
            )
            beta = -np.arcsin(-np.cos(cosines[0]))

        self.ends = [self.middle + L / 2 * nv, self.middle - L / 2 * nv]
        self.rotation = [alpha, beta, 0]


class DosePoint(object):
    """Point at which dose is calcluated.

    Attributes:
        x, y, z: x, y and z coordinates of point (cm).
        coords: array containing (x,y,z) coordinates (cm).
        name: name of dose point.
        ref: DICOM dose point reference number.
        dose: calculated dose (cGy).
        tpsdose: TPS dose obtained from DICOM file (cGy).
    """

    def __init__(self, coords, source, plan, name="", ref=""):
        """
        Args:
            coords: list of dose point coordinates
            rp: pydicom object of plan (RP) file.
            rs: pydicom object of structure (RS) file.
            name: name of dose point.
            ref: DICOM dose point reference number.
        """
        self.x, self.y, self.z = coords
        self.coords = np.array(coords)
        self.name = name
        self.ref = ref
        self.source = source
        if self.ref != "":
            self.get_tpsdose(plan.rp)
        self.calc_dose(source, plan)

    def __repr__(self):
        return self.name

    def get_tpsdose(self, rp):
        """Get dose calculated by TPS for this point.

        Args:
            rp: pydicom object of plan (RP) file.
        """
        tpsdose = 0

        if Tag(0x300A, 0x26) in rp[0x300A, 0x10][0].keys():
            for r in rp[0x300A, 0x10]:
                if self.ref == r[0x300A, 0x12].value:
                    tpsdose += r[0x300A, 0x26].value
        else:
            for cath in rp[0x300A, 0x230][0][0x300A, 0x280]:
                dwells = cath[0x300A, 0x2D0]

                if Tag(0x300C, 0x55) in dwells[-1].keys():
                    for r in dwells[-1][0x300C, 0x55]:
                        if self.ref == r[0x300C, 0x51].value:
                            tpsdose += r[0x300A, 0x10C].value

            tpsdose *= rp.FractionGroupSequence[0][0x300C, 0xA][0][0x300A, 0xA4].value
            tpsdose *= rp.FractionGroupSequence[0][0x300A, 0x78].value

            if rp.BrachyTreatmentType == "PDR":
                tpsdose *= rp[0x300A, 0x230][0][0x300A, 0x280][0][0x300A, 0x28A].value

        self.tpsdose = tpsdose

    def calc_dose(self, source, plan):
        """Calculate dose for this point

        Args:
            source: Source object.
            plan: Plan object.
        """
        dose = 0
        for d in plan.dwells:
            if d.t > 0:
                r = euclidzip(self.coords, d.middle)
                r1 = euclidzip(self.coords, d.ends[1])
                r2 = euclidzip(self.coords, d.ends[0])
                theta = np.arccos(
                    np.dot(
                        (d.ends[0] - d.ends[1]) / source.L, (self.coords - d.middle) / r
                    )
                )
                theta1 = np.arccos(
                    np.dot(
                        (d.ends[0] - d.ends[1]) / source.L,
                        (self.coords - d.ends[1]) / r1,
                    )
                )
                theta2 = np.arccos(
                    np.dot(
                        (d.ends[0] - d.ends[1]) / source.L,
                        (self.coords - d.ends[0]) / r2,
                    )
                )
                if (theta < 0.003) or ((np.pi - theta) < 0.003):
                    G = 1 / (r ** 2 - source.L ** 2 / 4)
                else:
                    B = np.abs(theta2 - theta1)
                    G = B / (source.L * r * np.sin(theta))

                dose += (
                    source.Sk
                    * source.Delta
                    * (G / source.G0)
                    * source.g(r)
                    * source.F(r, np.degrees(theta))
                    * d.t
                    / 3600
                )

        self.dose = dose / 100


class Plan(object):
    """Plan parameters.

    Attributes:
        rp: pydicom object of plan (RP) file.
        ROIs: list of ROI objects in plan.
        dwells: list of Dwell objects in plan.
    """

    def __init__(self, source, rp, rs, rd=None):
        """
        Args:
            source: Source object.
            rp: pydicom object of plan (RP) file.
            rs: pydicom object of structure (RS) file.
        """
        self.rp = rp
        if rp.BrachyTreatmentType == "PDR":
            self.frac = rp[0x300A, 0x230][0][0x300A, 0x280][0][0x300A, 0x28A].value
        elif rp.BrachyTreatmentType == "HDR":
            self.frac = rp.FractionGroupSequence[0][0x300A, 0x78].value
        self.get_ROIs(rs, rp, rd)
        self.get_dwells(source, rp)
        self.rx = rp.FractionGroupSequence[0][0x300C, 0xA][0][0x300A, 0xA4].value

    def get_ROIs(self, rs, rp, rd=None):
        """Get all structures in plan.

        Args:
            rs: pydicom object of structure (RS) file.
        """
        self.ROIs = []
        if rp.Manufacturer == "Nucletron":
            for roi in rp[0x300F, 0x1000][0][0x3006, 0x39]:
                if len(roi.ContourSequence) == 1:
                    self.ROIs.append(ROI(roi[0x3006, 0x84].value, None, rs, rp, rd))

        for struct in rs.StructureSetROISequence:
            self.ROIs.append(ROI(struct.ROINumber, struct.ROIName, rs, rp, rd))

    def get_dwells(self, source, rp):
        """Get all dwell points in plan.

        Args:
            source: Source object.
            rp: pydicom object of plan (RP) file.
        """
        c = 0
        self.dwells = []
        for cath in rp[0x300A, 0x230][0][0x300A, 0x280]:
            dwell_pts = cath[0x300A, 0x2D0]
            weight = cath[0x300A, 0x2C8].value
            total = cath[0x300A, 0x286].value * self.frac

            for roi in self.ROIs:
                if cath.ReferencedROINumber == roi.number:
                    app = roi
                elif (
                    rp.Manufacturer == "Nucletron"
                    and cath[0x300B, 0x1000].value == roi.number
                ):
                    app = roi

            for d in dwell_pts:
                x, y, z = d[0x300A, 0x2D4].value
                w = d[0x300A, 0x2D6].value - c
                if weight == 0:
                    t = 0
                else:
                    t = w / weight * total
                c = d[0x300A, 0x2D6].value
                self.dwells.append(Dwell([x / 10, y / 10, z / 10], t, source.L, app))


class ROI(object):
    """DICOM structure.

    Attributes:
        number: DICOM reference number.
        name: structure name.
        coords: array of (x,y,z) co-ordinates defining the structure (cm).
        tpsdvh: TPS-calculated cumulative DVH for this structure.
        dvh: cumulative DVH for this structure.
        dvhpts: co-ordinates for DVH calculation of this structure.
    """

    def __init__(self, number, name, rs, rp, rd=None):
        """
        Args:
            number: DICOM reference number.
            name: structure name.
            rs: pydicom object of structure (RS) file.
            rp: pydicom object of plan (RP) file.
        """
        self.number = number
        self.name = name
        self.coords = np.empty((0, 3))
        self.coordslist = []
        self.get_transform(rd)
        self.get_coords(rs, rp)

    def get_transform(self, rd):
        """Get transformation function to transform rotated orientation to default.

        Args:
            rd: pydicom object of dose (RD) file.
        """
        if rd:
            nnx = np.array(rd.ImageOrientationPatient[:3])
            nny = -np.array(rd.ImageOrientationPatient[3:])
            nnz = -np.cross(nnx, nny)
        else:
            nnx, nny, nnz = np.array([1, 0, 0, 0, 1, 0, 0, 0, 1], dtype=float).reshape(
                3, -1
            )

        nox, noy, noz = np.array([1, 0, 0, 0, 1, 0, 0, 0, 1], dtype=float).reshape(
            3, -1
        )

        top = [np.dot(nnx, n) for n in [nox, noy, noz]]
        mid = [np.dot(nny, n) for n in [nox, noy, noz]]
        bot = [np.dot(nnz, n) for n in [nox, noy, noz]]
        A = np.vstack((top, mid, bot))
        Ai = np.linalg.inv(A)

        def newcoords(coords):
            return np.vstack([np.dot(A, pt) for pt in coords])

        def oldcoords(coords):
            return np.vstack([np.dot(Ai, pt) for pt in coords])

        self.newcoords = newcoords
        self.oldcoords = oldcoords

    def get_coords(self, rs, rp):
        """Get co-ordinates for this structure from the RS file.

        Args:
            rs: pydicom object of structure (RS) file.
            rp: pydicom object of plan (RP) file.
        """

        if rs.Manufacturer == "Nucletron" and self.name == None:
            for roi in rp[0x300F, 0x1000][0][0x3006, 0x39]:
                if (
                    len(roi.ContourSequence) == 1
                    and roi[0x3006, 0x84].value == self.number
                ):
                    sli = roi.ContourSequence[0]
                    self.coords = np.append(
                        self.coords, np.array(sli.ContourData).reshape((-1, 3)), axis=0
                    )
        else:
            for roi in rs.ROIContourSequence:
                if (
                    roi[0x3006, 0x84].value == self.number
                    and Tag(0x3006, 0x40) in roi.keys()
                ):
                    for sli in roi.ContourSequence:
                        xyz = np.array(sli.ContourData).reshape((-1, 3))
                        xyz = np.round(self.newcoords(xyz), 1)
                        self.coordslist.append(xyz)
        if len(self.coordslist) > 0:
            self.coords = np.concatenate(self.coordslist)

    def get_TPS_DVH(self, rp, rs, rd):
        """Compute DVH for TPS-calculated dose distribution.

        Args:
            rp = pydicom object of plan (RP) file.
            rs = pydicom object of structure (RS) file.
            rd = pydicom object of dose (RD) file.
        """
        rx = rp.FractionGroupSequence[0][0x300C, 0xA][0][0x300A, 0xA4].value
        ix, iy, iz = self.newcoords([rd.ImagePositionPatient])[0]
        xcoords = ix + np.arange(0, rd.Columns * rd.PixelSpacing[0], rd.PixelSpacing[0])
        ycoords = iy - np.arange(0, rd.Rows * rd.PixelSpacing[1], rd.PixelSpacing[1])
        zcoords = np.round(iz + np.array(rd.GridFrameOffsetVector), 2)

        dose_ref = (rd.pixel_array * rd.DoseGridScaling)[:, :, :]

        dvh = {}
        coords = [(x, y) for y in ycoords for x in xcoords]

        for i in np.arange(len(rs.StructureSetROISequence)):
            ROIname = rs.StructureSetROISequence[i].ROIName
            if ROIname == self.name:
                contour = rs.ROIContourSequence[i]
                if "ContourSequence" in contour:
                    bool_ref = np.zeros(dose_ref.shape, dtype=bool)
                    dvh = []
                    for sli in contour.ContourSequence:
                        contourCoords = self.newcoords(
                            np.array(sli.ContourData).reshape((-1, 3))
                        )
                        index = np.round(np.mean(contourCoords[:, 2]), 2)
                        try:
                            boolslice = bool_ref[list(zcoords).index(index), :, :]
                            cPath = Path(contourCoords[:, :2])
                            inPath = cPath.contains_points(coords).reshape(
                                (dose_ref.shape[1], dose_ref.shape[2])
                            )
                            boolslice = np.logical_xor(inPath, boolslice)
                            bool_ref[list(zcoords).index(index), :, :] = boolslice
                        except:
                            pass
                    dvh = dose_ref[bool_ref == True]
                    self.tpsmin = dvh.min
                    self.tpsmax = dvh.max
                    self.tpsmean = dvh.mean
                    if len(dvh) > 0:
                        n, bins = np.histogram(dvh, 100, range=(0, rx * 10))
                        dvh = np.cumsum(n[::-1])[::-1]
                        dvh = dvh / dvh.max() * 100
                        self.tpsdvh = np.column_stack((bins[:-1], dvh))

    def get_DVH_pts(self, grid=2.5):
        """Calculate DVH calculation co-ordinates for this structure.
        Generates a grid of points at the specified resolution within
        the structure.

        Args:
            grid: grid size for calculation (default 2.5mm)
        """
        self.dvhpts = []
        slices = sorted(list(set(self.coords[:, 2])))
        for sli in slices:
            shapes = [x for x in self.coordslist if x[0, 2] == sli]
            if shapes:
                coords = np.concatenate(shapes)[:, :2]
                minx = coords[:, 0].min() - 0.5
                maxx = coords[:, 0].max() + 0.5
                miny = coords[:, 1].min() - 0.5
                maxy = coords[:, 1].max() + 0.5
                calcgrid = [
                    [x, y]
                    for y in np.arange(miny, maxy, grid)
                    for x in np.arange(minx, maxx, grid)
                ]
                boolgrid = [False] * len(calcgrid)
                for cont in shapes:
                    cPath = Path(cont[:, :2])
                    inPath = cPath.contains_points(calcgrid)
                    boolgrid = np.logical_xor(boolgrid, inPath)
                calcpts = [calcgrid[i] for i in np.where(boolgrid)[0]]
                self.dvhpts.extend(
                    [
                        list(self.oldcoords([[pt[0], pt[1], sli]])[0, [0, 1, 2]] / 10)
                        for pt in calcpts
                    ]
                )

    def __repr__(self):
        if self.name != None:
            return self.name
        else:
            return "Oncentra Applicator"


def euclidzip(v1, v2):
    """ Fast euclidean distance between two vectors v1 and v2. """
    dist = [(a - b) ** 2 for a, b in zip(v1, v2)]
    return math.sqrt(sum(dist))


def fastinterp(xx, yy):
    """ 1D interpolation. """

    def interpout(x):
        I = np.searchsorted(xx, x)
        x1 = xx[I - 1]
        x2 = xx[I]

        y1 = yy[I - 1]
        y2 = yy[I]

        return y1 + (y2 - y1) * (x - x1) / (x2 - x1)

    return interpout


from bisect import bisect_left


def bilinearinterp(xi, yi, values):
    """ Bilinear interpolation. """

    def interpolate(x, y):
        # local lookups
        i = bisect_left(xi, x) - 1
        j = bisect_left(yi, y) - 1

        x1, x2 = xi[i : i + 2]
        y1, y2 = yi[j : j + 2]
        z11, z12 = values[j][i : i + 2]
        z21, z22 = values[j + 1][i : i + 2]

        return (
            z11 * (x2 - x) * (y2 - y)
            + z21 * (x - x1) * (y2 - y)
            + z12 * (x2 - x) * (y - y1)
            + z22 * (x - x1) * (y - y1)
        ) / ((x2 - x1) * (y2 - y1))

    return interpolate
