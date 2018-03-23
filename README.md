# pyTG43

pyTG43 uses user-defined source specification data as well as TPS-defined dwell positions and times to calculate dose at a specific point using the 2D formalism outlined in [AAPM TG-43](http://dx.doi.org/10.1118/1.1646040).

pyTG43 has only been tested on files from Eclipse and Oncentra treatment planning systems. No warranty or support is provided.

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

All source specification data can be downloaded in spreadsheet form from the [ESTRO website](https://www.estro.org/about/governance-organisation/committees-activities/tg43). Just specify the directory in which you have placed the HDR and/or PDR spreadsheets, and pyTG43 will extract the relevant data.

## Example calculations

Examples are proivded for both Eclipse and Oncentra treatment planning sytems, and can be run using `python examples.py`

```bash
$ python examples.py
----------------------------------------
 Eclipse (HDR)
----------------------------------------
┌───────────┬───────┬──────┬───────┬───────────┬────────────┬──────────┐
│ Name      │ X     │ Y    │ Z     │ TPS (cGy) │ Calc (cGy) │ Diff (%) │
├───────────┼───────┼──────┼───────┼───────────┼────────────┼──────────┤
│ PtA_left  │ 1.91  │ 2.28 │ -1.25 │ 600.156   │ 600.302    │ -0.024   │
│ PtA_right │ -2.09 │ 2.11 │ -1.25 │ 613.616   │ 613.782    │ -0.027   │
└───────────┴───────┴──────┴───────┴───────────┴────────────┴──────────┘
----------------------------------------
 Oncentra (PDR)
----------------------------------------
┌──────┬───────┬────────┬───────┬───────────┬────────────┬──────────┐
│ Name │ X     │ Y      │ Z     │ TPS (cGy) │ Calc (cGy) │ Diff (%) │
├──────┼───────┼────────┼───────┼───────────┼────────────┼──────────┤
│ P1   │ -1.04 │ -62.60 │ -3.47 │ 440.772   │ 440.499    │ 0.062    │
│ P2   │ -0.23 │ -62.60 │ -3.43 │ 438.617   │ 438.459    │ 0.036    │
│ -    │ -2.27 │ -63.70 │ -2.42 │ 173.064   │ 173.054    │ 0.006    │
│ -    │ 1.70  │ -63.74 │ -2.60 │ 127.385   │ 127.605    │ -0.173   │
│ -    │ -3.66 │ -58.42 │ -5.19 │ 25.354    │ 25.681     │ -1.291   │
│ -    │ -0.06 │ -58.42 │ -5.19 │ 37.546    │ 37.727     │ -0.483   │
│ -    │ -3.57 │ -59.25 │ -4.72 │ 33.656    │ 33.896     │ -0.714   │
│ -    │ -3.29 │ -60.13 │ -4.34 │ 43.770    │ 43.910     │ -0.319   │
│ -    │ -3.04 │ -61.03 │ -4.00 │ 52.495    │ 52.455     │ 0.077    │
│ -    │ -2.80 │ -61.95 │ -3.68 │ 65.908    │ 65.679     │ 0.348    │
│ -    │ 0.43  │ -59.25 │ -4.72 │ 42.213    │ 42.318     │ -0.248   │
│ -    │ 0.71  │ -60.13 │ -4.34 │ 48.794    │ 48.802     │ -0.017   │
│ -    │ 0.96  │ -61.03 │ -4.00 │ 57.933    │ 57.815     │ 0.203    │
│ -    │ 1.20  │ -61.95 │ -3.68 │ 73.067    │ 72.862     │ 0.280    │
└──────┴───────┴────────┴───────┴───────────┴────────────┴──────────┘
```
