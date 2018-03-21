import pyTG43
import pydicom

rp = pydicom.dcmread('example/RP.example.dcm')
rs = pydicom.dcmread('example/RS.example.dcm')

pyTG43.tpsComp(rp, rs, 'example/')