# HETR-Flight-Dynamics-LQR
Nonlinear 6-DOF flight dynamics simulation and LQR stability augmentation for a Hybrid Electric Tiltrotor (XV-15 derivative)

---

## Overview

This repository contains the complete MATLAB flight dynamics and control
Implementation for an XV-15 derivative Hybrid Electric Tiltrotor (HETR)


The analysis covers:
- Nonlinear 6-DOF rigid-body simulation (1000-second cruise phase)
- Analytical trim computation across the conversion corridor (110–190 m/s)
- Numerical linearisation and open-loop modal analysis (5 natural modes)
- LQR Stability Augmentation System design (longitudinal + lateral)
- Gain-scheduled LQR across the full conversion corridor (214–369 KCAS)
- Nonlinear closed-loop validation with disturbance rejection simulation

---

## Project Highlights

- Developed a nonlinear 6-DOF rigid-body flight-dynamics model
- Performed trim and stability analysis across the conversion corridor
- Designed longitudinal and lateral LQR stability augmentation
- Implemented gain scheduling across 213.8–369.4 KCAS
- Evaluated open-loop and closed-loop modal characteristics
- Validated nonlinear disturbance rejection
---

## Tools

- MATLAB
- Control System Toolbox
- Numerical ODE integration
- State-space modelling
- LQR control design
---

## Aircraft Parameters

| Parameter | Value |
|-----------|-------|
| Mass (m) | 7,800 kg |
| Wing area (S) | 15.433 m² |
| Wingspan (b) | 9.8 m |
| Max thrust (T_max) | 8,000 N |
| Cruise altitude | 15,000 ft (ρ = 0.6564 kg/m³) |
| Cruise speed (Va,trim) | 157.4 m/s (306.0 KCAS) |
| Trim L/D | 16.42 |

---

## Key Results

### Open-Loop Modal Analysis

| Mode | ωn (rad/s) | ζ | Level |
|------|-----------|---|-------|
| Short-period | 0.584 | 0.509 | Level 1 ✓ |
| Phugoid | 0.090 | 0.120 | Level 1 ✓ |
| Dutch roll | 2.065 | 0.028 | Level 2 (needs SAS) |
| Roll subsidence | 0.349 | 1.000 | Level 2 |
| Spiral | 0.059 | 1.000 | Level 1 ✓ |

### LQR Closed-Loop Results

| Mode | OL ζ | CL ζ | Improvement |
|------|------|------|-------------|
| Short-period | 0.509 | 0.726 | +43% |
| Phugoid | 0.120 | 0.705 | +488% |
| Dutch roll | 0.028 | 0.680 | ×24 |

All 5 closed-loop modes achieve **ADS-33E Level 1** handling qualities.

### Disturbance Rejection
- Perturbation applied: +5° pitch, +3° roll simultaneously
- Full state recovery to trim: **< 8 seconds**
- All control surfaces within **±25° saturation limit**

### Gain Schedule
Closed-loop short-period ζ maintained in band **0.718–0.728** (Level 1)
across the full conversion corridor from 213.8 to 369.4 KCAS.

---

## Simulation Results

### Open-Loop Cruise Simulation
![Open-loop simulation]
<img width="1170" height="723" alt="open_loop_simulation" src="https://github.com/user-attachments/assets/67b49a14-8218-4623-a4c2-ba6407a2e492" />


### Open-Loop vs LQR Closed-Loop Disturbance Rejection
![OL vs CL comparison]
<img width="1920" height="926" alt="open_vs_closed_loop" src="https://github.com/user-attachments/assets/faabc19c-42c9-4ef1-b075-a22f06aa7fe5" />


### LQR Control Surface Deflections
![Control surfaces]<img width="1920" height="926" alt="LQR_control_surface_deflections" src="https://github.com/user-attachments/assets/9b5189e1-7c74-4be4-b53f-7a2b6d8990a3" />


### Eigenvalue Migration — Open-Loop to Closed-Loop
![Eigenvalue map]<img width="1920" height="926" alt="Combine_eihen_plot" src="https://github.com/user-attachments/assets/68e81782-7021-4230-80d6-68c03518028b" />


---



## References

- NASA SP-2000-4517 — XV-15 Tilt Rotor Research Aircraft
- Nelson, R.C. — *Flight Stability and Automatic Control* (2nd ed.)
- Stevens, Lewis & Johnson — *Aircraft Simulation and Control* (3rd ed.)
- ADS-33E-PRF — Handling Qualities Requirements for Military Rotorcraft
