# PrPOD-ROM FEMU for Online Dynamic Analysis

MATLAB implementation of a **frequency-domain finite element model updating (FEMU)** framework that integrates a  
**progressive proper orthogonal decomposition (PrPOD)-based reduced-order model (ROM)** with **particle swarm optimization (PSO)**.  
By reducing repeated full-model simulations, the framework enables efficient model updating and fast online dynamic response analysis.  

---

## 📌 Related Publication

This repository provides MATLAB codes for the following paper:

**Chanwoo Lee**, Namsu Jeon, and Hyung-Jo Jung  
*Rapid finite element model updating for online dynamic analysis via progressive frequency proper orthogonal decomposition*  
Engineering Structures (2026)  
DOI: [https://doi.org/10.1016/j.engstruct.2026.122228](https://www.sciencedirect.com/science/article/pii/S0141029626001410)

---

## 🎯 Motivation

FEMU typically requires iterative simulations, and high-fidelity FE models can be computationally prohibitive for online use.  
This work proposes an **offline–online** strategy:

- **Offline phase**: construct a **parametric frequency POD-ROM** using FRF data.
- **Online phase**: perform **FEMU + online dynamic analysis** efficiently using the ROM.

The method is designed to improve feasibility for rapid decision-making after events such as earthquakes.

<p align="left">
  <img src="figures/overview.png" width="700">
</p>

---

## ⚙️ Algorithm Overview

Core components:

- **Frequency-domain FEMU** using FRF-based comparison (robust without explicit mode matching)
- **Frequency POD-ROM** for projection-based reduced governing equations
- **PrPOD (Progressive POD)**: incrementally increases POD mode count during optimization
- **PSO (Particle Swarm Optimization)** with adaptive inertia and multi-criteria convergence checks

The workflow follows an **offline–online** structure.

<p align="left">
  <img src="figures/framework.png" width="700">
</p>

### Workflow

- **Input** → Operating points {μ^(i)}, system matrices M(μ), C(μ), K(μ), frequency grid, measured FRFs, PSO settings, and progressive mode counts (r1 < r2 < ... < M).
- **Offline (ROM construction)**  
  1) Compute FRFs at operating points  
  2) Build snapshot matrix from selected frequency samples  
  3) Truncated SVD → POD bases  
  4) Congruence / Procrustes alignment for consistent interpolation  
- **Online (FEMU + dynamic analysis)**  
  1) Interpolate aligned bases for trial parameter μ_trial  
  2) Assemble reduced operators via Galerkin projection  
  3) Compute reduced FRF and reconstruct full response at measured DOFs  
  4) Evaluate FRF-based objective and update parameters via PSO  
  5) Progressively increase rank r until final rank M  
  6) With updated μ*, perform fast time-domain dynamic analysis using ROM  

<p align="left">
  <img src="figures/workflow.png" width="520">
</p>

---

## 🚀 Running the Code

### 1) Offline: Build parametric ROM
Run:
- `main/main_offline_build_ROM.m`  
  → Computes FRFs at operating points, builds frequency POD bases, and performs basis alignment.

### 2) Online: Model updating (PSO + PrPOD)
Run:
- `main/main_online_FEMU_PSO.m`  
  → Executes progressive PSO-based FEMU with increasing POD ranks.

### 3) Online: Fast time-domain response (after updating)
Run:
- `main/main_online_time_response.m`  
  → Computes online dynamic response using the final updated ROM.

---

## 📌 Case Studies

### (A) 5-Story Frame Structure
Validated with a five-story frame model (numerical + experimental).  
Demonstrates accurate parameter estimation and response reproduction even under challenging conditions (e.g., closely spaced modes, sensor noise, sparse measurements).

### (B) Nuclear Containment Building (APR1400)
Demonstrated on a high-fidelity nuclear containment FE model under seismic excitation to show scalability and online applicability.

---

## 📄 Citation

If you use this repository in your research, please cite the related paper above.
