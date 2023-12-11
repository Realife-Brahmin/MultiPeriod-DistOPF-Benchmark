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

## Procedure of Creating these .DSS files
A MATLAB optimization program is used to solve for Optimal Power Flow (OPF) on the system 
with some user-specified inputs, such as degree of GED penetration, load profile, 
PV profile, etc.
After the OPF run, all decision variables, along with system details, are
input into a separate script which interfaces with the OpenDSS MATLAB engine,
to solve for powerflow, to validate MATLAB's results.
The model once compiled, is then 'saved' as a Master.dss file and other .dss
files required by it to run the powerflow.

It may be noted that GEDs such as PVs and Batteries are modelled simply as
generators with a fixed real/reactive power value equal to the corresponding 
decision variable output by the OPF solver.

## How to Run Powerflow for the Model using OpenDSS GUI:

Just open Master-OpenDSS.dss (or Master.dss if you wanna specify your own solver
instructions) and Ctrl+A Ctrl+D.