# Multi-Period Distributed Optimal Power Flow
## Optimization for Balanced Three-Phase Power Distribution Networks with Renewables and Storage in MATLAB. 
Naive Brute Force Multi-Period OPF. A spatially decomposed, temporally brute-forced MPOPF has been implemented.

Objectves currently covered: 
- Loss Minimization

## Description of the Modelling of the Radial Power Distribution System

### Description of State Variables

| Variable Notation | Variable Description                  | Number of Variables | Nature of Constraint |
|-------------------|---------------------------------------|---------------------|----------------------|
| $P^{t}_{ij}$      | Real Power flowing in branch          | $m$                 | Nonlinear            |
| $Q^{t}_{ij}$      | Reactive Power flowing in branch      | $m$                 | Nonlinear            |
| $l^{t}_{ij}$      | Square of Magnitude of branch Current | $m$                 | Nonlinear            |
| $v^{t}_{j}$       | Square of Magnitude of node Voltage   | $N$                 | Nonlinear            |
| $B^{t}_{j}$       | Battery State of Charge               | $n_{B}$             | Linear               |

### Description of Control Variables

| Variable Notation | Variable Description                     | Number of Variables | Nature of Constraint |
|-------------------|------------------------------------------|---------------------|----------------------|
| $q^{t}_{D_j}$     | Reactive Power of DER (via inverter)     | $n_{D}$             | Linear*              |
| $P^{t}_{c_j}$     | Charging Power of Battery                | $n_{B}$             | Linear               |
| $P^{t}_{d_j}$     | Discharging Power of Battery             | $n_{B}$             | Linear               |
| $q^{t}_{B_j}$     | Reactive Power of Battery (via inverter) | $n_{B}$             | Linear*              |

### Description of Independent Variables

| Variable Notation | Variable Description            | Number of Variables | Nature of Constraint |
|-------------------|---------------------------------|---------------------|----------------------|
| $P^{t}_{L_j}$     | Real Power Demand               | $N$                 | Linear               |
| $Q^{t}_{L_j}$     | Reactive Power Demand           | $N$                 | Linear               |
| $P^{t}_{D_j}$     | Real Power of DER               | $n_{D}$             | Linear*              |
| $B^{0}_{j}$       | Battery Initial State of Charge | $n_{B}$             | Linear               |


