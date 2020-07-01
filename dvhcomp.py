import pydicom
from itertools import cycle
from pyTG43 import *
from bokeh.plotting import figure, show, output_file
from bokeh.palettes import Category20
from bokeh.models import CrosshairTool, HoverTool, Legend

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

cmap = Category20[20]
itr = cycle(cmap)

output_file("dvh.html")

p = figure(plot_width=900, plot_height=560, x_range=(0,plan.rx*10), y_range=(0,101), active_scroll='wheel_zoom', toolbar_location='above')
p.xaxis.axis_label = 'Dose (Gy)'
p.yaxis.axis_label = 'Relative Volume (%)'

items = []
for roi in plan.ROIs:
    if hasattr(roi, 'dvh') and hasattr(roi, 'tpsdvh'):
        clr = next(itr)
        a = p.line(roi.tpsdvh[:,0],roi.tpsdvh[:,1], line_color=clr, line_width=1.5)
        b = p.line(roi.dvh[:,0],roi.dvh[:,1], line_color=clr, line_dash=(2,2), line_width=1.5)
        items.append((roi.name, [a,b]))

p.toolbar.autohide = True
hover = HoverTool(
    tooltips = [
    ("Volume","$y{1.11}%"),
    ("Dose", "$x{1.11} Gy")
    ],
    show_arrow = False,
    line_policy = 'interp'
)
cross = CrosshairTool(
    line_alpha = 0.5
)

legend = Legend(
    items = items,
    location=(10,0),
    click_policy='hide',
    background_fill_alpha=1
)
p.add_layout(legend, 'right')

legend2 = Legend(
    items = [
        ('pyTG43', [p.line([0],[0],line_color='black',line_dash=(2,2))]),
        ('TPS', [p.line([0],[0],line_color='black')])
    ],
    location='top_right',
    background_fill_alpha=0,
    border_line_alpha=0
)
p.add_layout(legend2, 'center')

p.add_tools(cross, hover)

show(p)