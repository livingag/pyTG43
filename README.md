# pyTG43

pyTG43 uses user-defined source specification data as well as TPS-defined dwell positions and times to calculate dose at a specific point using the 2D formalism outlined in [AAPM TG-43](http://dx.doi.org/10.1118/1.1646040).

pyTG43 has only been tested on files from Varian brachytherapy systems and is distributed with no warranty or support included.

## Installation

```
git clone https://github.com/livingag/pyTG43.git
cd pyTG43
python setup.py install
```

## Calculation details

The co-ordinate system used by the software is [IEC 61217](ftp://dicom.nema.org/MEDICAL/dicom/2015b/output/chtml/part03/figures/PS3.3_C.8.8.25-2a.svg), which is used by most planning systems. Keep in mind this is different from the [DICOM](http://dicom.nema.org/DICOM/2013/output/chtml/part17/figures/PS3.17_FFF.1.2-3.svg) co-ordinate system.

Dwell point co-ordinates and times are extracted from the DICOM plan file, along with the applicator structures from the DICOM structure set file.

For each dwell point, the orientation of the line source is determined by the vector of the two nearest points of the applicator structure for that dwell point. The relevant angles and distances required by the TG-43 formalism are then easily determined through vector calculus.

## Source specification data

Data specifying the source dose rate constant, anisotropy factor, and radial dose function is required. An example set of data (not corresponding to any real source) is shown in the `example/` directory. Only one of each source type can be input at any time.

A `.yaml` file is used to specify all required source data. The example shows the two available methods to specifying the radial dose function; with an n-th order polynomial, or with a table of input data to be interpolated.

```
delta: dose rate constant
anisotropy:
  filename: filename for data containing anisotropy data
radial:
  coeffs: coefficients of polynomial fit (list)
  filename: filename for data containing radial data
```
An error indicating that source length data is missing means that the source length is not present in the DICOM plan file (the case for some TPSes). Specify the source length (in cm) the `sourcespec.yaml` file using the `length:` attribute.

## Example calculation

Example taken from Varian BrachyVision system with tandem and ovoid applicator in a square block of water. Can be run using `python example.py`

```python
>> import pyTG43
>> import pydicom

>> rp = pydicom.dcmread('example/RP.example.dcm')
>> rs = pydicom.dcmread('example/RS.example.dcm')

>> pyTG43.tpsComp(rp, rs, 'example/')

┌───────────┬───────┬──────┬───────┬───────────┬────────────┬──────────┐
│ Name      │ X     │ Y    │ Z     │ TPS (cGy) │ Calc (cGy) │ Diff (%) │
├───────────┼───────┼──────┼───────┼───────────┼────────────┼──────────┤
│ PtA_left  │ 1.91  │ 2.28 │ -1.25 │ 600.156   │ 602.332    │ -0.363   │
│ PtA_right │ -2.09 │ 2.11 │ -1.25 │ 613.616   │ 615.860    │ -0.366   │
└───────────┴───────┴──────┴───────┴───────────┴────────────┴──────────┘
```
