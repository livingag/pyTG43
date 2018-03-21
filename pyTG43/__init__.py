import numpy as np
import pydicom
from pydicom.tag import Tag
from scipy.interpolate import interp1d, interp2d
from scipy.spatial.distance import euclidean
import yaml
from terminaltables import SingleTable

def tpsComp(rp, rs, directory):

    points = []

    for p in rp[0x300a, 0x10]:
        if Tag(0x300a,0x18) in p.keys():
            x, z, y = p[0x300a, 0x18].value
            name = p[0x300a, 0x16].value
            ref = p[0x300a, 0x12].value
            points.append(DosePoint([x/10, y/10, z/10], rp, rs, name, ref, directory))

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


class Dwell(object):
    def __init__(self, coords, t, L, app):
        self.x, self.y, self.z = coords
        self.middle = np.array(coords)
        self.t = t
        self.get_source_ends(app, L)

    def get_source_ends(self, app, L):
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
    def __init__(self, coords, rp, rs, name='', ref='', directory='.'):
        self.x, self.y, self.z = coords
        self.coords = np.array(coords)
        self.name = name
        self.ref = ref
        self.get_tpsdose(rp)
        self.calc_dose(rp, rs, directory)

    def __repr__(self):
        return self.name

    def get_tpsdose(self, rp):
        tpsdose = 0

        if Tag(0x300a, 0x26) in rp[0x300a, 0x10][0].keys():
            for r in rp[0x300a, 0x10]:
                if self.ref == r[0x300a, 0x12].value:
                    tpsdose += r[0x300a, 0x26]
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

    def calc_dose(self, rp, rs, directory):
        r0 = 1
        theta0 = np.pi/2
        Sk = rp.SourceSequence[0].ReferenceAirKermaRate
        L = rp.SourceSequence[0].ActiveSourceLength/10

        with open(directory+'/sourcespec.yaml','r') as stream:
            sourcespec = yaml.load(stream)

        source = sourcespec[rp.BrachyTreatmentType]

        Delta = source['delta']

        Fi = np.loadtxt(directory+'/'+source['anisotropy']['filename'], delimiter=',')
        Fi = interp2d(Fi[0,:][1:],Fi[:,0][1:],Fi[1:,1:])
        def F(r,theta):
            if r > 15:
                return Fi(15,theta)
            else:
                return Fi(r,theta)

        if 'filename' in source['radial'].keys():
            g = np.loadtxt(directory+'/'+source['radial']['filename'], delimiter=',')
            g = interp1d(g[:,0],g[:,1])
        elif 'coeff' in source['radial'].keys():
            coeff = source['radial']['coeff']
            def g(r):
                if r < 0.15:
                    return g(0.15)
                elif r > 15:
                    return g(12)*np.exp((r-12)/(15-12)*(np.log(g(15))-np.log(g(12))))
                nth = len(coeff)
                out = 0
                for n in range(nth):
                    out += coeff[n]*r**n
                return out

        G0 = 2*np.arctan((L/2)/r0)/(L*r0*np.sin(theta0))

        ROIs = []
        for struct in rs.StructureSetROISequence:
            ROIs.append(ROI(struct.ROINumber, struct.ROIName, rs))

        time = 0
        dwells = []
        for cath in rp[0x300a, 0x230][0][0x300a, 0x280]:
            dwell_pts = cath[0x300a, 0x2d0]

            for roi in ROIs:
                if cath.SourceApplicatorID == roi.name:
                    app = roi

            for d in dwell_pts:
                x, z, y = d[0x300a, 0x2d4].value
                t = d[0x300a, 0x2d6].value - time
                time = d[0x300a, 0x2d6].value
                dwells.append(Dwell([x/10, y/10, z/10], t, L, app))

        dose = 0
        for d in dwells:
            if d.t > 0:
                r = euclidean(self.coords,d.middle)
                r1 = euclidean(self.coords,d.ends[1])
                r2 = euclidean(self.coords,d.ends[0])
                theta = np.arccos(np.dot((d.ends[0]-d.ends[1])/L,(self.coords - d.middle)/r))
                theta1 = np.arccos(np.dot((d.ends[0]-d.ends[1])/L,(self.coords - d.ends[1])/r1))
                theta2 = np.arccos(np.dot((d.ends[0]-d.ends[1])/L,(self.coords - d.ends[0])/r2))
                if (theta < 0.003) or ((np.pi-theta) < 0.003):
                    G = 1/(r**2 - L**2/4)
                else:
                    B = np.abs(theta2 - theta1)
                    G = B/(L*r*np.sin(theta))

                dose += Sk * Delta * (G/G0) * g(r) * F(r,np.degrees(theta))[0] * d.t / 3600

        self.dose = dose


class ROI(object):
    def __init__(self, number, name, rs):
        self.number = number
        self.name = name
        self.coords = np.empty((0,3))
        self.get_coords(rs)

    def get_coords(self, rs):
        for roi in rs.ROIContourSequence:
            if roi[0x3006,0x84].value == self.number and Tag(0x3006, 0x40) in roi.keys():
                for sli in roi.ContourSequence:
                    self.coords = np.append(self.coords,np.array(sli.ContourData).reshape((-1, 3)),axis=0)
                self.coords /= 10
                self.coords = self.coords[:, [0, 2, 1]]

    def __repr__(self):
        return self.name
