# Multi-Period Distributed Optimal Power Flow
## Optimization for Balanced Three-Phase Power Distribution Networks with Renewables and Storage in MATLAB. 
Naive Brute Force Multi-Period OPF. A spatially decomposed, temporally brute-forced MPOPF has been implemented.

Objectves currently covered: 
- Loss Minimization

## Description of the Modelling of the Radial Power Distribution System

### Description of State Variables

| Variable Notation | Variable Description                                      | Number of Variables | Nature of Constraint |
|-------------------|-----------------------------------------------------------|---------------------|----------------------|
| $P_{ij}$          | Real Power flowing in branch $(i, j) \in \mathbb{L}$      | $m$                 | Nonlinear            |
| $Q_{ij}$          | Reactive Power flowing in branch $(i, j) \in \mathbb{L}$  | $m$                 | Nonlinear            |
| $l_{ij}$          | Square of Magnitude of Current flowing in branch $(i, j)$ | $m$                 | Nonlinear            |
| $v_{j}$           | Square of Magnitude of Voltage at node $j \in \mathbb{N}$ | $N$                 | Nonlinear            |
| $B_{j}$           | Battery State of Charge at node $j \in \mathbb{B}$        | $n_{B}$             | Linear               |

### Description of Control Variables

| Variable Notation | Variable Description                                                | Number of Variables | Nature of Constraint |
|-------------------|---------------------------------------------------------------------|---------------------|----------------------|
| $q_{D_j}$         | Reactive Power of DER (via inverter) at node $j \in \mathbb{D}$     | $n_{D}$             | Linear*              |
| $P_{c_j}$         | Charging Power of Battery at node $j \in \mathbb{B}$                | $n_{B}$             | Linear               |
| $P_{d_j}$         | Discharging Power of Battery at node $j \in \mathbb{B}$             | $n_{B}$             | Linear               |
| $q_{B_j}$         | Reactive Power of Battery (via inverter) at node $j \in \mathbb{B}$ | $n_{B}$             | Linear*              |

### Description of Independent Variables

| Variable Notation | Variable Description                                       | Number of Variables | Nature of Constraint |
|-------------------|------------------------------------------------------------|---------------------|----------------------|
| $P_{L_j}$         | Real Power Demand at node $j \in \mathbb{N}$               | $N$                 | Linear               |
| $Q_{L_j}$         | Reactive Power Demand at node $j \in \mathbb{N}$           | $N$                 | Linear               |
| $P_{D_j}$         | Real Power of DER at node $j \in \mathbb{D}$               | $n_{D}$             | Linear*              |
| $B^{0}_{j}$       | Battery Initial State of Charge at node $j \in \mathbb{B}$ | $n_{B}$             | Linear               |


