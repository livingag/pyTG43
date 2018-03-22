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

The co-ordinate system used by the software is [IEC 61217](https://i.imgur.com/k926EqO.png), which is used by most planning systems. Keep in mind this is different from the [DICOM](http://dicom.nema.org/DICOM/2013/output/chtml/part17/figures/PS3.17_FFF.1.2-3.svg) co-ordinate system.

Dwell point co-ordinates and times are extracted from the DICOM plan file, along with the applicator structures from the DICOM structure set file.

For each dwell point, the orientation of the line source is determined by the vector of the two nearest points of the applicator structure for that dwell point. The relevant angles and distances required by the TG-43 formalism are then easily determined through vector calculus.

## Source specification data

All source specification data can be downloaded in spreadsheet form from the [ESTRO website](https://www.estro.org/about/governance-organisation/committees-activities/tg43). Just specify the directory in which you have placed the HDR and/or PDR spreadsheets.

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
