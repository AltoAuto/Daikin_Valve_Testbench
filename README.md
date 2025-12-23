<h1 align = 'center'> Daikin Valve Test Bench Digital Twin

<p align="center">
  <img src="assets/poster_thumb.png" width="950" alt="Design Poster Overview">
</p>

## What this project is 
This repository contains a **Simscape-based digital twin** of a vapor-compression refrigerant loop designed to support **expansion valve testing** (EXV/TXV). The model is built to answer the same questions an HVAC OEM would ask:

- Can the bench hit the required **ṁ–Δp valve test window** without exceeding pressure limits?
- Does **superheat control** have authority across operating points (no saturation/windup)?
- Does the system remain within **safe envelopes** under worst-case conditions?
- Are **water-side loads** and **energy balance** physically consistent?
- Does the loop behave safely during **startup and load transients**?

---
## System architecture

- **Main model:** `model/MainBench.slx`
- **Control:** EXV PI superheat regulation (and TXV behavior where applicable)
- **Validation style:** phased “test suite” (component → open-loop → closed-loop → safety → energy → dynamics)

---
## Validation framework

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
