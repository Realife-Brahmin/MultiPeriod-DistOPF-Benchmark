# System Description

- System Name: IEEE 123 Bus System

- Network Model: Balanced Three-Phase (Single-Phase)

- Description: This IEEE 123 Bus System is a Radial Balanced Three-Phase distribution system with a single substation. Some nodes also have (optionally) reactive-controlled photovoltaics and reactive-controlled batteries on them.

## How to Run Powerflow for the Model using OpenDSS GUI:

Just open Master-OpenDSS.dss in OpenDSS and Ctrl+A Ctrl+D.

## High level overview of the repository
ieee123_1ph/
│
├── archives/                 # Archive files, such as previous versions and backups
│
├── components/               # Modules for different electrical components
│   ├── PV_configs/           # Configuration files for Photovoltaic systems
│   │   ├── PV10.dss          # Configuration for PV system 10
│   │   ├── PV30.dss          # Configuration for PV system 30
│   │   ├── PV50.dss          # Configuration for PV system 50
│   │   ├── PV100.dss         # Configuration for PV system 100
│   │   └── PV100_reactive.dss # Configuration for reactive components of PV system 100
│   ├── Battery_configs/      # Configuration files for Battery systems
│   │   ├── Battery10.dss     # Configuration for Battery 10
│   │   └── StorageControl10.dss # Configuration for the storage controller of Battery 10
│   ├── Capacitor.dss           # Capacitor modules
│   └── LineCode.dss  # contains definitions of some line modules used for defining branchData in IEEE123SinglePh.dss
│   ├── Regulator.dss           # Voltage regulator modules
│   ├── RegulatorControl.dss    # Regulator Controller modules IEEE123SinglePh.DSS
│   └── IEEE123SinglePh.dss     # busData and branchData
│   └── [Other Component Files] # Additional component configuration files
│
├── data/                     # Input data for the simulation
│   └── LoadShape1.CSV    # CSV file containing load shape data
│   └── [Other Data Files]    # Additional data files
│
├── extra/                    # Extra resources and miscellaneous files
│
├── resources/                # Documentation and reference materials
│
├── results/                  # Output data from simulations, including monitor data and plots
├── pv010_batt010/            # Example configuration folder for a specific PV and Battery setup
│   ├── counts*.csv           # CSV files containing counts data
│   ├── loads*.csv            # CSV files containing loads data
│   ├── powers*.csv           # CSV files containing power data
│   ├── summary*.csv          # CSV files containing summary data
│   └── voltages*.csv         # CSV files containing voltages data
│
├── scripts/                  # Custom scripts that are called by the Master-OpenDSS.dss file
│   └── solveAndExportFor24h.dss # Script for solving and exporting results over 24 hours
│
│
└── Master-OpenDSS.dss        # Main OpenDSS model file for the ieee123_1ph project
