## System Description

- System Name: IEEE 123 Bus System

- Network Model: Balanced Three-Phase (Single-Phase)

- Description: This IEEE 123 Bus System is a Radial Balanced Three-Phase distribution system with a single substation. For this test case, some nodes also have reactive-controlled photovoltaics on them. Batteries are currently NOT modelled in it, but will be added in a later version.

## Components Modelled/Yet-to-be-Modelled: 

### Note on Components Modelled: 
All Grid Edge Devices (GEDs) such as photovoltaics (PVs) and Batteries are modelled as Generators. These generators can have fixed negative real/reactive powers. The assigned powers are as dictated by the optimization routine.

### Components Already Modelled:
- PV (photovoltaics), with two-quadrant reactive control

### Components Yet to be Modelled: 
- Batteries with two-quadrant reactive power control

## Current Contents part of the IEEE123 Balanced Three-Phase (Single-phase) OpenDSS model:

- BusCoords.dss (Automatically created by OpenDSS)
- BusVoltageBases.dss (Setting the Base Voltage value for the system)
- Generator.dss (to model all PVs)
- GrowthShape.dss (Automatically created by OpenDSS)
- Line.dss (Branch Data)
- Load.dss (Contains all Power Demands)
- LoadShape.dss (Automatically created by OpenDSS)
- Master.dss
- Master-OpenDSS.dss (Same as Master.dss, but with a few added lines, can just be Ctrl+A Ctrl+D to get powerflow results)
- Spectrum.dss (Automatically created by OpenDSS)
- TCC_Curve.dss (Automatically created by OpenDSS)
- VSource.dss (Substation node data)

## How to Run Powerflow for the Model using OpenDSS GUI:

Just open Master-OpenDSS.dss (or Master.dss if you wanna specify your own solver instructions) and Ctrl+A Ctrl+D.