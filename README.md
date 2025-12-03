# **Daikin HVAC Predictive Refrigeration Loop — Digital Twin (Simscape) **

*A geometry-anchored, EXV/TXV-validated, component-correct digital twin for an R-32 refrigeration loop.*

---

## **1. Project Overview**

This repository contains a **validated**, **multi-phase**, **predictive** HVAC refrigeration-loop model built in **MATLAB Simscape Fluids**.
The goal is to demonstrate:

* Predictive behavior from geometry (not empirical “black-box” UA)
* Accurate EXV/TXV mass-flow + superheat control
* Full-cycle refrigerant energetics
* Component-level validation
* A test bench system usable by future students/industry

This repo serves as:
* A **modeling reference**
* A **validation package**
* A **starting point** for future predictive HVAC research

---

## **2. Repository Structure**

```
/ComponentTestHarness/
    - Compressor harness (map-based or PD-based)
    - TL-2P Evaporator/Condenser harness
    - EXV/TXV flow harness
    - Pipe/accumulator models
    - Validation plots & test scripts

/MainModel/
    - Full-cycle refrigerant model
    - Geometry-anchored heat exchangers
    - EXV/TXV control layers
    - System-level scripts
    - Operating point solver
    - Transient tests
```

Each folder represents a **phase of validation**.

---

## **3. Model Architecture**

High-level block diagram:

* Compressor
* Condenser (TL-2P predictive model)
* EXV
* Evaporator (TL-2P predictive model)
* Controllers (superheat control, compressor speed, etc.)

All major elements are built with:

* **Simscape Two-Phase Fluids**
* **Thermal Liquid ↔ 2-Phase predictive exchanger** (see reference:
  *Condenser Evaporator (TL-2P)* model )

This enables:
* Zone-based L/M/V tracking
* Fin geometry
* Pressure drop modeling
* Predictive UA (no black box)

---

## **4. Component Validation Summary**

The repo includes data + scripts for:

### Compressor

* Map ingestion (Daikin-provided CSV)
* Speed sweep validation
* Mass flow vs pressure ratio consistency

### EXV/TXV

* Flow coefficient mapping
* Superheat response
* Compressor-sweep SH mapping
* Linearity + hysteresis tests

### Heat Exchangers (Predictive TL-2P)

* Geometry-based UA
* E-NTU
* Zone fraction evolution
* ΔP modeling vs expected correlations

Reference model description:
(See MathWorks heat exchanger detailed documentation) 

---

## **5. Full System Model**

System-level simulation features:

* Evaporator/condenser geometry
* R-32 property functions
* Charge inventory tracking
* Variable compressor speed
* Disturbance tests
* Load sweeps
* Fouling scenarios

---

## **6. Validation Workflow **

This README should list the same **phased validation pipeline** we’re building:

1. **Phase 0 — Requirements & Acceptance Criteria**
2. **Phase 1 — Component Harness Validation**
3. **Phase 2 — Open-Loop Cycle Validation**
4. **Phase 3 — Closed-Loop EXV/TXV SH Control**
5. **Phase 4 — Safety Envelopes & Fault Cases**
6. **Phase 5 — Water Side & AHRI 550/590 compliance**
7. **Phase 6 — Dynamic startup/load-step validation**
8. **Phase 7 — Final Validation Report & Traceability**

Each phase has its own folder in `/plots`, `/data`, or `/scripts` later.

---

## **7. How to Run the Model**

1. Install MATLAB R2025b+ with:

   * Simscape Fluids
   * Simscape Electrical
   * Simscape Foundation Library

2. Open:

```
/MainModel/RefrigerantLoop_main.slx
```

3. Run:

```
run startup.m
```

4. Then:

```
sim('RefrigerantLoop_main')
```

---

## 📎 **8. Included Documents**

You can list uploaded engineering artifact PDFs here:

* System schematic diagram 
* Predictive TL-2P heat exchanger reference 
* Daikin Applied Expansion Valve specifications
* Test stand design report
* Compressor map documentation
* AHRI 550/590 reference standard for cycle performance 

This builds credibility.

---

## **9. FAQ / Known Limitations**

Examples:

* Current compressor map is from public Daikin R-32 dataset
* TXV model validated only under steady superheat
* Oil return not yet included
* Charge model is quasi-static

---

## **10. How to Cite / Use This Repo**

“Please cite this repository if you use it for research or teaching.”
