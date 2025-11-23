# Chaboche Model Calibration Using Particle Swarm Optimization (MATLAB)

This repository provides MATLAB code for calibrating a **nonlinear Chaboche plasticity model** to experimental stress–strain data.  
The workflow implements:

- Multi-component kinematic hardening (X₁, X₂, X₃)  
- Isotropic hardening (R)  
- A **β-switch rule** for the third backstress component  
- Elastic modulus estimation from initial experimental data  
- Parameter identification using **Particle Swarm Optimization (PSO)**  
- Model vs. experiment comparison plots  

The script reads experimental CSV data, optimizes 9 material parameters, and simulates stress response over the full loading path.

---

## Features

### ✔️ Reads experimental data  
`1.5percent.csv` containing:  
- Total strain  
- Plastic strain  
- Stress (MPa)

### ✔️ Automatic estimation of Young’s modulus  
Using a linear fit to initial elastic region.

### ✔️ Chaboche hardening model
Includes:


- **Kinematic hardening (3 components):**  
  $\( X_i' = \frac{2}{3} C_i \, dp \, n - \gamma_i X_i \, dp \)$

- **Isotropic hardening:**  
  $\( R' = b(Q - R) dp \)$

### ✔️ PSO-based parameter identification  
Optimizes:
[C1, g1, C2, g2, C3, g31, g32, Q, b]

### ✔️ Visualization  
Plots model prediction vs. experimental curve.

---

## Usage

### 1. Place your CSV file  
Expected format (3 columns):


### 2. Run the script in MATLAB

```matlab
run your_script_name.m



