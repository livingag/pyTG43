import pyTG43
import pydicom

rp = pydicom.dcmread('examples/eclipse/RP.eclipse.dcm')
rs = pydicom.dcmread('examples/eclipse/RS.eclipse.dcm')

print('-'*40)
print(' Eclipse (HDR)')
print('-'*40)
pyTG43.tpsComp(rp, rs, 'examples/eclipse/')

rp = pydicom.dcmread('examples/oncentra/RP.oncentra.dcm')
rs = pydicom.dcmread('examples/oncentra/RS.oncentra.dcm')

print('-'*40)
print(' Oncentra (PDR)')
print('-'*40)
pyTG43.tpsComp(rp, rs, 'examples/oncentra/')