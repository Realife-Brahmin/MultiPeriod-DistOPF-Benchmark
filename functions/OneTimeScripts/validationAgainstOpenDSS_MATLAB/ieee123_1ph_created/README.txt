# System Description

- System Name: IEEE 123 Bus System

- Network Model: Balanced Three-Phase (Single-Phase)

- Description: This IEEE 123 Bus System is a Radial Balanced Three-Phase distribution system with a single substation. Some nodes also have (optionally) reactive-controlled photovoltaics and reactive-controlled batteries on them.

## How to Run Powerflow for the Model using OpenDSS GUI:

Just open Master.dss in OpenDSS and Ctrl+A Ctrl+D.

## High level overview of the repository
ieee123_1ph_OpenDSS
│
├── components
│   ├── Battery_configs
│   │   ├── Battery10.dss
│   │   ├── StorageControl10.dss
│   │   └── ...
│   │
│   ├── Monitors
│   │   └── voltageMonitors10.dss
│   │
│   ├── PVconfigs
│   │   ├── PV10.dss
│   │   ├── PV30.dss
│   │   └── ...
│   │
│   ├── Capacitor.dss
│   ├── IEEE123SinglePh.DSS
│   ├── LineCode.DSS
│   ├── Regulator.dss
│   └── RegulatorControl.dss
│
├── data
│   ├── LoadShape_PV.xlsx
│   ├── LoadShape_Storage.xlsx
│   └── ...
│
├── openDSSFilesCreation
│   └── ...
│
├── results
│   ├── figures
│   │   └── ...
│   │
│   ├── pv000_batt000
│   └── ...
│
├── scripts
│   ├── plotLoadShape.jl
│   ├── plotMaster.jl
│   ├── solveAndExportForOneHour.dss
│   └── ...
│
├── Master.dss
├── OpenDSS-files-for-MPOPF-verification.code-workspace
├── README.md
└── README.txt

