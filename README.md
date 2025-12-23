<h1 align="center">Daikin Valve Test Bench Digital Twin
<hr>
<p align="center">
  <img src="assets/poster_thumb.png" width="950" alt="Design Poster Overview">
</p>

## What This Project Is 
This repository contains a **Simscape-based digital twin** of a vapor-compression refrigerant loop designed to support **expansion valve testing** (EXV/TXV). The model is built to answer the same questions an HVAC OEM would ask:

- Can the bench hit the required **ṁ–Δp valve test window** without exceeding pressure limits?
- Does **superheat control** have authority across operating points (no saturation/windup)?
- Does the system remain within **safe envelopes** under worst-case conditions?
- Are **water-side loads** and **energy balance** physically consistent?
- Does the loop behave safely during **startup and load transients**?

---
## System Architecture

- **Main model:** `model/MainBench.slx`
- **Control:** EXV PI superheat regulation (and TXV behavior where applicable)
- **Validation style:** phased “test suite” (component → open-loop → closed-loop → safety → energy → dynamics)

---
## Validation Framework

The validation is organized into sequential phases with pass/fail criteria and documented evidence.

- **Phase 0 — Requirements & acceptance criteria**
- **Phase 1 — Component harness validation**
- **Phase 2 — Open-loop cycle validation**
- **Phase 3 — Closed-loop SH control validation**
- **Phase 4 — Safety & worst-case envelope**
- **Phase 5 — Water-side sizing & energy balance**
- **Phase 6 — Dynamic transient validation (startup / load step)**
- **Phase 8 — Final validation report + traceability**
---

## Repo Structure
```
daikin-valve-testbench/
├─ docs/
│  ├─ Validation_Manual.pdf
│  ├─ Design_Poster.pdf
├─ README.md
├─ startup.m
│  └─ Schematic.png
├─ model/
│  └─ MainBench.slx
├─ cases/
│  └─ case_nominal.m
├─ scripts/
│  └─ run_main.m
├─ assets/
│  └─ poster_thumb.png
└─ results/        (gitignored)
```

---
## Authorship & Attribution

This repository is maintained by **Aidan W.** as a portfolio-quality technical package derived from a senior design project context.

**Primary contributions**
- System architecture and model decomposition
- Simscape-based component and loop modeling
- Validation methodology and acceptance criteria
- Control integration and test-case execution
- Technical documentation and presentation assets

**Team context**
- University senior design collaboration  
- Contributor names are listed in the design poster and final report where applicable
