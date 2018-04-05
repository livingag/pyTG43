from glob import glob

import numpy as np
import pydicom
from pydicom.tag import Tag
from scipy.interpolate import interp1d, interp2d
from scipy.spatial.distance import euclidean
from terminaltables import SingleTable
from matplotlib.path import Path
import xlrd

def tpsComp(rp, rs, directory):
    """Calculate and compare dose at reference points with TPS.

    Args:
        rp: pydicom object of RP file.
        rs: pydicom object of RS file.
        directory: directory containing source data.

    Returns:
        points: list of DosePoint objects.
    """

    points = []

    source = Source(rp, directory)
    plan = Plan(source, rp, rs)

    for p in rp[0x300a, 0x10]:
        if Tag(0x300a,0x18) in p.keys():
            x, z, y = p[0x300a, 0x18].value
            name = p[0x300a, 0x16].value
            ref = p[0x300a, 0x12].value
            points.append(DosePoint([x/10, y/10, z/10], source, plan, name, ref, directory))

    table_data = [['Name','X','Y','Z','TPS (cGy)','Calc (cGy)','Diff (%)']]

    for pt in points:
        table_data.append([pt.name,
                           '{:.2f}'.format(pt.x),
                           '{:.2f}'.format(pt.y),
                           '{:.2f}'.format(pt.z),
                           '{:.3f}'.format(pt.tpsdose),
                           '{:.3f}'.format(pt.dose),
                           '{:.3f}'.format((1-(pt.dose/pt.tpsdose))*100)])

    print(SingleTable(table_data).table)
    return points


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
            rp: pydicom object of RP file.
            directory: directory containing source data.
        """
        r0 = 1
        theta0 = np.pi/2

        self.Sk = rp.SourceSequence[0].ReferenceAirKermaRate

        fname = glob(directory+'/*'+rp.BrachyTreatmentType.lower()+'*.xls')[0]
        wb = xlrd.open_workbook(fname)
        sh = wb.sheets()[-1]

        self.L = sh.row(9)[2].value
        self.Delta = sh.row(4)[2].value

        col = 5
        while sh.row(10)[col].ctype == 2:
            col += 1

        Fi = np.ones((1, col-5))
        gi = np.ones((1, 2))
        Ft = []
        Fr = [r.value for r in sh.row(10)[5:col]]

        for row in np.arange(11,sh.nrows):
            if sh.row(row)[1].ctype != 0:
                gi = np.vstack([gi,np.array([r.value for r in sh.row(row)[1:3]])])
            Fi = np.vstack([Fi,np.array([r.value for r in sh.row(row)[5:col]])])
            Ft.append(sh.row(row)[4].value)

        Fi = Fi[1:,:]
        gi = gi[1:,:]

        self.Fi = interp2d(Fr,Ft,Fi)
        self.gi = interp1d(gi[:,0],gi[:,1])

        self.G0 = 2*np.arctan((self.L/2)/r0)/(self.L*r0*np.sin(theta0))

    def F(self, r, theta):
        """Anisotropy function.

        Args:
            r: distance (cm).
            theta: angle (radians).
        """
        if r > 10:
            return self.Fi(10,theta)
        else:
            return self.Fi(r,theta)
    def g(self, r):
        """Radial dose function.

        Args:
            r: distance (cm).
        """
        if r < 0.15:
            return self.gi(0.15)
        elif r > 10:
            return self.gi(8)*np.exp((r-8)/(10-8)*(np.log(self.gi(10))-np.log(self.gi(8))))
        else:
            return self.gi(r)


class Dwell(object):
    """Planned dwell point.

    Attributes:
        x, y, z: x, y and z coordinates (cm).
        middle: (x,y,z) coordinates of middle of source (cm).
        t: dwell time (s).
        ends: list containing coordinates of each endd of the source (cm).
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
        self.get_source_ends(app, L)

    def get_source_ends(self, app, L):
        """Calculates coordinates of the ends of the source.

        Args:
            app: applicator ROI object dwell point is a part of.
            L: source length (cm).
        """
        closestd = [+np.inf,+np.inf]
        closest = [0,0]
        for pt in app.coords:
            d = euclidean(self.middle,pt)
            if d < closestd[0]:
                closestd[0] = d
                closest[0] = pt
            elif d < closestd[1]:
                closestd[1] = d
                closest[1] = pt

        v = closest[0] - closest[1]
        r = euclidean(closest[0],closest[1])
        nv = v / r

        self.ends = [self.middle + L/2 * nv, self.middle - L/2 * nv]


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
    def __init__(self, coords, source, plan, name='', ref=''):
        """
        Args:
            coords: list of dose point coordinates
            rp: pydicom object of RP file.
            rs: pydicom object of RS file.
            name: name of dose point.
            ref: DICOM dose point reference number.
        """
        self.x, self.y, self.z = coords
        self.coords = np.array(coords)
        self.name = name
        self.ref = ref
        self.source = source
        if self.ref != '':
            self.get_tpsdose(plan.rp)
        self.calc_dose(source, plan)

    def __repr__(self):
        return self.name

    def get_tpsdose(self, rp):
        """Get dose calculated by TPS for this point.

        Args:
            rp: pydicom object of RP file.
        """
        tpsdose = 0

        if Tag(0x300a, 0x26) in rp[0x300a, 0x10][0].keys():
            for r in rp[0x300a, 0x10]:
                if self.ref == r[0x300a, 0x12].value:
                    tpsdose += r[0x300a, 0x26].value
        else:
            for cath in rp[0x300a, 0x230][0][0x300a, 0x280]:
                dwells = cath[0x300a, 0x2d0]

                if Tag(0x300c, 0x55) in dwells[-1].keys():
                    for r in dwells[-1][0x300c, 0x55]:
                        if self.ref == r[0x300c, 0x51].value:
                            tpsdose += r[0x300a, 0x10c].value

            tpsdose *= rp.FractionGroupSequence[0][0x300c, 0xa][0][0x300a, 0xa4].value

            if rp.BrachyTreatmentType == 'PDR':
                tpsdose *= rp[0x300a, 0x230][0][0x300a, 0x280][0][0x300a, 0x28a].value

        self.tpsdose = tpsdose * 100

    def calc_dose(self, source, plan):
        """Calculate dose for this point

        Args:
            source: Source object.
            plan: Plan object.
        """
        dose = 0
        for d in plan.dwells:
            if d.t > 0:
                r = euclidean(self.coords,d.middle)
                r1 = euclidean(self.coords,d.ends[1])
                r2 = euclidean(self.coords,d.ends[0])
                theta = np.arccos(np.dot((d.ends[0]-d.ends[1])/source.L,(self.coords - d.middle)/r))
                theta1 = np.arccos(np.dot((d.ends[0]-d.ends[1])/source.L,(self.coords - d.ends[1])/r1))
                theta2 = np.arccos(np.dot((d.ends[0]-d.ends[1])/source.L,(self.coords - d.ends[0])/r2))
                if (theta < 0.003) or ((np.pi-theta) < 0.003):
                    G = 1/(r**2 - source.L**2/4)
                else:
                    B = np.abs(theta2 - theta1)
                    G = B/(source.L*r*np.sin(theta))

                dose += source.Sk * source.Delta * (G/source.G0) * source.g(r) * source.F(r,np.degrees(theta))[0] * d.t / 3600

        self.dose = dose

class Plan(object):
    """Plan parameters.

    Attributes:
        rp: pydicom object of RP file.
        ROIs: list of ROI objects in plan.
        dwells: list of Dwell objects in plan.
    """
    def __init__(self, source, rp, rs):
        """
        Args:
            source: Source object.
            rp: pydicom object of RP file.
            rs: pydicom object of RS file.
        """
        self.rp = rp
        self.get_ROIs(rs, rp)
        self.get_dwells(source, rp)
        self.rx = rp.FractionGroupSequence[0][0x300c, 0xa][0][0x300a, 0xa4].value

        if rp.BrachyTreatmentType == 'PDR':
            self.rx *= rp[0x300a, 0x230][0][0x300a, 0x280][0][0x300a, 0x28a].value
    def get_ROIs(self, rs, rp):
        """Get all structures in plan.

        Args:
            rs: pydicom object of RS file.
        """
        self.ROIs = []
        if rp.Manufacturer == 'Nucletron':
            for roi in rp[0x300f, 0x1000][0][0x3006, 0x39]:
                if len(roi.ContourSequence) == 1:
                    self.ROIs.append(ROI(roi[0x3006,0x84].value,None,rs,rp))

        for struct in rs.StructureSetROISequence:
            self.ROIs.append(ROI(struct.ROINumber, struct.ROIName, rs, rp))
    def get_dwells(self, source, rp):
        """Get all dwell points in plan.

        Args:
            source: Source object.
            rp: pydicom object of RP file.
        """
        c = 0
        self.dwells = []
        for cath in rp[0x300a, 0x230][0][0x300a, 0x280]:
            dwell_pts = cath[0x300a, 0x2d0]
            weight = cath[0x300a, 0x2c8].value
            total = cath[0x300a, 0x286].value

            for roi in self.ROIs:
                if cath.SourceApplicatorID == roi.name:
                    app = roi
                elif rp.Manufacturer == 'Nucletron' and \
                     cath[0x300b,0x1000].value == roi.number:
                    app = roi

            for d in dwell_pts:
                x, z, y = d[0x300a, 0x2d4].value
                w = d[0x300a, 0x2d6].value - c
                if weight == 0:
                    t = 0
                else:
                    t = w/weight * total
                c = d[0x300a, 0x2d6].value
                self.dwells.append(Dwell([x/10, y/10, z/10], t, source.L, app))


class ROI(object):
    """DICOM structure.

    Attributes:
        number: DICOM reference number.
        name: structure name.
        coords: array of (x,y,z) co-ordinates defining the structure (cm).
        dvh: cumulative DVH for this structure.
    """
    def __init__(self, number, name, rs, rp):
        """
        Args:
            number: DICOM reference number.
            name: structure name.
            rs: pydicom object of RS file.
            rp: pydicom object of RP file.
        """
        self.number = number
        self.name = name
        self.coords = np.empty((0,3))
        self.get_coords(rs,rp)

    def get_coords(self, rs, rp):
        """Get co-ordinates for this structure from the RS file.

        Args:
            rs: pydicom object of RS file.
            rp: pydicom object of RP file.
        """

        if rs.Manufacturer == 'Nucletron' and self.name == None:
            for roi in rp[0x300f, 0x1000][0][0x3006, 0x39]:
                if len(roi.ContourSequence) == 1 and roi[0x3006, 0x84].value == self.number:
                    sli = roi.ContourSequence[0]
                    self.coords = np.append(self.coords,np.array(sli.ContourData).reshape((-1, 3)),axis=0)
        else:
            for roi in rs.ROIContourSequence:
                if roi[0x3006,0x84].value == self.number and Tag(0x3006, 0x40) in roi.keys():
                    for sli in roi.ContourSequence:
                        self.coords = np.append(self.coords,np.array(sli.ContourData).reshape((-1, 3)),axis=0)

        self.coords /= 10
        self.coords = self.coords[:, [0, 2, 1]]

    def get_DVH(self,source,plan,grid=0.2):
        """Calculate cumulative DVH for this structure.

        Args:
            source: source object.
            plan: plan object.
            grid: grid size for calculation (default 2mm)
        """
        dvh = []
        if len(self.coords) > 50:
            slices = sorted(list(set(self.coords[:,1])))
            sthick = slices[1] - slices[0]
            for sli in slices:
                xy = self.coords[np.where(self.coords[:,1] == sli)][:,(0,2)]
                minx = xy[:,0].min()
                maxx = xy[:,0].max()
                miny = xy[:,1].min()
                maxy = xy[:,1].max()
                calcgrid = [[x, y] for y in np.arange(miny,maxy,grid) for x in np.arange(minx,maxx,grid)]
                cPath = Path(xy)
                inPath = cPath.contains_points(calcgrid)
                calcpts = [calcgrid[i] for i in np.where(inPath)[0]]
                for pt in calcpts:
                    pt = [pt[0], sli, pt[1]]
                    dvh.append(DosePoint(pt,source,plan).dose)

            self.volume = len(dvh) * (grid * grid * sthick)
            n, bins = np.histogram(dvh,100,range=(0,plan.rx*11))
            dvh = np.cumsum(n[::-1])[::-1]
            dvh = dvh / dvh.max() * 100
            self.dvh = np.column_stack((bins[:-1],dvh))
    def __repr__(self):
        if self.name != None:
            return self.name
        else:
            return 'Oncentra Applicator'
