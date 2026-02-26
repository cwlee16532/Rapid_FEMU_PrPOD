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
DOI: https://doi.org/10.1016/j.engstruct.2026.122228  

> The published Version of Record is available at the journal website via the DOI link above.

---

## 🎯 Motivation

Finite element model updating (FEMU) typically requires iterative simulations, and high-fidelity FE models can be computationally prohibitive for online use.  

This work proposes an **offline–online strategy**:

- **Offline phase** → Construct a parametric frequency POD-ROM using FRF data  
- **Online phase** → Perform FEMU and fast dynamic analysis efficiently using the ROM  

The objective is to significantly reduce computational cost while maintaining accuracy in parameter estimation and response prediction.

---

## ⚙️ Algorithm Overview

Core components of the framework include:

- **Frequency-domain FEMU** using FRF-based comparison (robust without explicit mode matching)  
- **Frequency POD-ROM** for projection-based reduced governing equations  
- **PrPOD (Progressive POD)** → Incrementally increases POD mode count during optimization  
- **PSO (Particle Swarm Optimization)** with adaptive inertia and multi-criteria convergence checks  

The algorithm follows an **offline–online structure**.

### Workflow

**Input**  
Operating points {μ^(i)}, system matrices M(μ), C(μ), K(μ), frequency grid, measured FRFs, PSO settings, and progressive mode counts (r1 < r2 < ... < M).

---

### Offline Phase (ROM Construction)

1. Compute FRFs at predefined operating points  
2. Build snapshot matrix from selected frequency samples  
3. Perform truncated SVD to extract POD bases  
4. Apply congruence / Procrustes alignment for consistent basis interpolation  

---

### Online Phase (FEMU + Dynamic Analysis)

1. Interpolate aligned bases for trial parameter μ_trial  
2. Assemble reduced operators via Galerkin projection  
3. Compute reduced FRF and reconstruct full response at measured DOFs  
4. Evaluate FRF-based objective function and update parameters via PSO  
5. Progressively increase POD rank r until final rank M  
6. Perform fast time-domain dynamic analysis using the updated ROM  

---

## 🚀 Running the Code

### 1) Offline: Build Parametric ROM

Run:
- `main/main_offline_build_ROM.m`  

This script:
- Computes FRFs at operating points  
- Builds frequency POD bases  
- Performs basis alignment for interpolation  

---

### 2) Online: Model Updating (PSO + PrPOD)

Run:
- `main/main_online_FEMU_PSO.m`  

This script:
- Executes progressive PSO-based FEMU  
- Increases POD rank progressively  
- Identifies updated structural parameters  

---

### 3) Online: Fast Time-Domain Response

Run:
- `main/main_online_time_response.m`  

This script:
- Performs online dynamic response analysis  
- Uses the final updated reduced-order model  

---

## 📌 Case Study: 5-Story Frame Structure

The framework is validated using a laboratory-scale five-story frame structure.

The case study demonstrates:

- Accurate parameter identification under stiffness variations  
- Robust FRF matching without explicit mode pairing  
- Efficient computation through progressive ROM expansion  
- Reliable time-domain response prediction after updating  

The results confirm that the PrPOD-ROM FEMU framework significantly reduces computational cost while maintaining high accuracy.


---

## ⚖️ Copyright Notice

This repository contains original MATLAB implementation codes developed by the author.  
The published journal article is not redistributed here.  
