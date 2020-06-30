import pydicom
from itertools import cycle
from pyTG43 import *
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

rp = pydicom.dcmread('examples/PDR/RP.PDR.dcm')
rs = pydicom.dcmread('examples/PDR/RS.PDR.dcm')
rd = pydicom.dcmread('examples/PDR/RD.PDR.dcm')

source = Source(rp, 'examples/PDR/')
plan = Plan(source, rp, rs, rd)

for roi in plan.ROIs:
    if roi.name in ['ctv','bladder','rectum']:
        roi.get_DVH_pts()
        roi.get_TPS_DVH(rp,rs,rd)

calcDVHs(source,plan,plan.rx*10,['ctv','bladder','rectum'])

plt.rcParams['figure.figsize'] = [10, 7]
cmap = plt.get_cmap('tab20').colors
itr = cycle(cmap)

i = 0
plt.gca().set_ylim([0,101])
plt.gca().set_xlim([0,plan.rx*10])
for roi in plan.ROIs:
    if hasattr(roi, 'dvh'):
        clr = next(itr)
        plt.plot(roi.tpsdvh[:,0],roi.tpsdvh[:,1],c=clr,linestyle='-')
        plt.plot(roi.dvh[:,0],roi.dvh[:,1],c=clr,linestyle='--')

plt.xlabel('Dose (Gy)')
plt.ylabel('Relative Volume (%)')
ax = plt.gca().add_artist(plt.legend(plt.gca().lines[::2],[roi.name for roi in plan.ROIs if hasattr(roi, 'dvh')],bbox_to_anchor=(1.04,1), loc="upper left", borderaxespad=0))
solid_line = mlines.Line2D([], [], color='black', linestyle='-',label='TPS')
dashed_line = mlines.Line2D([], [], color='black', linestyle='--',label='pyTG43')
plt.legend(handles = [solid_line,dashed_line],bbox_to_anchor=(1.04,0), loc="lower left", borderaxespad=0)