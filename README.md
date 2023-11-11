# Multi-Period Distributed Optimal Power Flow
## Optimization for Balanced Three-Phase Power Distribution Networks with Renewables and Storage in MATLAB. 
Naive Brute Force Multi-Period OPF. A spatially decomposed, temporally brute-forced MPOPF has been implemented.

Objectves currently covered: 
- Loss Minimization

## Description of the Modelling

### Description of State Variables

| Variable Notation | Variable Description                                      | Number of Variables | Nature of Constraint |
|-------------------|-----------------------------------------------------------|---------------------|----------------------|
| $P_{ij}$          | Real Power flowing in branch $(i, j) \in \mathbb{L}$      | $m$                 | Nonlinear            |
| $Q_{ij}$          | Reactive Power flowing in branch $(i, j) \in \mathbb{L}$  | $m$                 | Nonlinear            |
| $l_{ij}$          | Square of Magnitude of Current flowing in branch $(i, j)$ | $m$                 | Nonlinear            |
| $v_{j}$           | Square of Magnitude of Voltage at node $j \in \mathbb{N}$ | $N$                 | Nonlinear            |
| $B_{j}$           | Battery State of Charge at node $j \in \mathbb{B}$        | $n_{B}$             | Linear               |

