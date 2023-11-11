# Multi-Period Distributed Optimal Power Flow
## Optimization for Balanced Three-Phase Power Distribution Networks with Renewables and Storage in MATLAB. 
Naive Brute Force Multi-Period OPF. A spatially decomposed, temporally brute-forced MPOPF has been implemented.

## Objectives currently covered: 
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
| $q^{t}_{D_j}$     | Reactive Power of DER (via inverter)     | $n_{D}$             | Linear<sup id="a1">[1](#f1)</sup>              |
| $P^{t}_{c_j}$     | Charging Power of Battery                | $n_{B}$             | Linear               |
| $P^{t}_{d_j}$     | Discharging Power of Battery             | $n_{B}$             | Linear               |
| $q^{t}_{B_j}$     | Reactive Power of Battery (via inverter) | $n_{B}$             | Linear<sup id="a1">[1](#f1)</sup>              |

### Description of Independent Variables

| Variable Notation | Variable Description            | Number of Variables | Nature of Constraint |
|-------------------|---------------------------------|---------------------|----------------------|
| $P^{t}_{L_j}$     | Real Power Demand               | $N$                 | Linear               |
| $Q^{t}_{L_j}$     | Reactive Power Demand           | $N$                 | Linear               |
| $P^{t}_{D_j}$     | Real Power of DER               | $n_{D}$             | Linear<sup id="a1">[1](#f1)</sup>              |
| $B^{0}_{j}$       | Battery Initial State of Charge | $n_{B}$             | Linear               |

### Miscellaneous Notation

| Variable Notation | Variable Description                                                               | Cardinality |
|-------------------|------------------------------------------------------------------------------------|-------------|
| $\mathbb{N}$      | Set of all the nodes                                                               | $N$         |
| $\mathbb{L}$      | Set containing all the branches                                                    | $m$         |
| $\mathbb{D}$      | Set containing all the nodes containing DERs. $\mathbb{D} \subset \mathbb{N}$      | $n_{D}$     |
| $\mathbb{B}$      | Set containing all the nodes containing Batteries. $\mathbb{B} \subset \mathbb{N}$ | $n_{B}$     |
| $\mathbb{T}$      | Set containing all the time-periods                                                | $T$         |
| $j$               | Denotes a node. $j \in \mathbb{N}$                                                 |             |
| $(i, j)$          | Denotes a branch connecting nodes $i$ and $j$. $(i, j) \in \mathbb{L}$             |             |
| $t$               | Denotes a time-period<sup id="a2">[2](#f2)</sup>. $t \in \mathbb{T}$                                        |             |

### Notes
1. <span id="f1"> Current modelling. Future modelling will incorporate reactive power as a non-linear function wrt maximum apparent power and real power. [↩](#a1)</span>
2. <span id="f2"> Except when used as a superscript in denoting Battery SOC $B^{t}_j$, $t$ refers to the average value of the variable _within_ the time-period $t$. For Battery SOC, $B^{t}_j$ refers to the value of SOC _at the end_ of time-period $t$. [↩](#a2)</span>
