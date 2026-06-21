# hydro-rmt-geysers-floods

**From Geysers to Megafloods: Random Matrix Theory Reveals Universal Level Repulsion in Hydrogeological Charge-and-Release Systems**

A single geyser exceeds GUE rigidity while multi-source mixing collapses to Poisson.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)

**Author:** Ruqing Chen · GUT Geoservice Inc., Montreal · ruqing@hotmail.com

---

## Core Finding

The **pressure-shadow hypothesis**: a fluid-release event (geyser eruption,
glacial-lake outburst) depletes local fluid pressure and thermal energy, which
must re-accumulate before the next event — imposing level repulsion analogous
to eigenvalue repulsion. Tested across three hydrogeological configurations
spanning **six orders of magnitude in timescale**:

| Target | Data source | Timescale | n | ⟨r⟩ | CV | Class |
|---|---|---|---|---|---|---|
| **A: Old Faithful Geyser** | Real (Azzalini 1990) | ~70 min | 271 | **0.738** | 0.20 | **Super-GUE** |
| B: Yellowstone Basin mixed | Forward model (15 sources) | ~0.6 min | 800 | 0.383 | 1.07 | Poisson |
| **C: N. Atlantic IRD events** | Real (DSDP 609, HSG peaks) | ~1.5 kyr | 34 | **0.550** | 0.77 | **GOE** |

Reference values: Poisson = 0.386 · GOE = 0.536 · GUE = 0.603

## The Super-GUE Regime — Key Result

Old Faithful's ⟨r⟩ = 0.738 is **the strongest repulsion observed in any
geological system** in the entire RMT program. It exceeds the GUE limit
(0.603) with Brody β = 3.0, reflecting the extreme confinement of a single
narrow hydrothermal conduit with a near-deterministic two-state charge–release
cycle. The geyser is *more rigid than random matrix eigenvalues* — approaching
the crystalline limit of perfect periodicity (β → ∞, CV → 0).

## Spectral Superposition — Largest Contrast

The multi-geyser basin (⟨r⟩ = 0.383) and Old Faithful (⟨r⟩ = 0.738) differ
by Δ⟨r⟩ = 0.355 — **the largest single-vs-mixed contrast** in the program.
Mixing 15 independent charge–release oscillators erases every trace of
single-source memory, exactly as the theorem predicts.

## Glacial Outburst Floods — GOE at Kiloyear Timescales

35 hematite-stained grain (HSG) peaks from DSDP Site 609 (Bond & Obrochta
2012) record ice-rafted debris pulses spanning 14.4–70.8 ka BP. These
glacial charge–release events (ice-sheet buildup → surge → depletion → rebuild)
show GOE-class repulsion (⟨r⟩ = 0.550), the same universality class as mantle
plumes and ore deposits operating at megayear timescales.

## Grand Synthesis — Eight Systems, One Principle

| System | Timescale | ⟨r⟩ | Class | Study |
|---|---|---|---|---|
| Old Faithful geyser | ~70 min | 0.738 | Super-GUE | **this work** |
| Tethyan porphyry Cu | ~1 Myr | 0.712 | GUE | Paper 4 |
| Orogenic gold | ~10 Myr | 0.678 | GUE | Paper 4 |
| Mantle plumes (×4) | ~5 Myr | 0.630 | GOE–GUE | Paper 3 |
| Andean porphyry Cu | ~2 Myr | 0.601 | GOE–GUE | Paper 4 |
| Nanling W–Sn | ~5 Myr | 0.574 | GOE | Paper 4 |
| N. Atl. IRD events | ~1.5 kyr | 0.550 | GOE | **this work** |

**Eight physically unrelated Earth systems, one principle: single long-memory
charge–release sources repel; superposed sources randomize.**

## Repository Structure
```
hydro-rmt-geysers-floods/
├── README.md · LICENSE · requirements.txt · CITATION.cff · .zenodo.json
├── paper/
│   ├── paper.tex · paper.pdf      # 12 pp.
│   └── figs/                      # figures embedded by LaTeX
├── code/
│   └── hydro_rmt_pipeline.py      # Three-target analysis + visualization
├── data/
│   └── 94-609_Site_HSG_IG_MIS4-2.tab   # DSDP 609 HSG (Bond/Obrochta 2012)
├── figures/                       # standalone PNGs (300 dpi)
│   ├── fig1_hydro_panel.png       # 1×3 spacing distributions
│   └── fig2_grand_synthesis.png   # 8-system bar chart
└── results/                       # JSON outputs
```

## Method
1. **Old Faithful**: 272 consecutive eruption waiting times (Azzalini & Bowman 1990)
   → normalize → compute ⟨r⟩, CV, Brody β, KS tests
2. **Basin mixed**: forward model of 15 independent Poisson geysers superposed
   → spectral superposition theorem validation
3. **IRD events**: DSDP 609 HSG record → peak detection (≥15%, prominence ≥4)
   → 35 peaks, 34 spacings in kiloyears → RMT analysis

## Reproduce
```bash
pip install -r requirements.txt
cd code
python hydro_rmt_pipeline.py
```

## Five-Racetrack RMT Program
1. Geological boundaries (Myr) → GOE — [zenodo 20766310](https://zenodo.org/records/20766310)
2. Seismotectonic rhythms → scale-dependent — [zenodo 20768130](https://zenodo.org/records/20768130)
3. Mantle plumes (Gyr) → single-source GOE — [zenodo 20768420](https://zenodo.org/records/20768420)
4. Metallogeny (Myr) → single ore system GOE/GUE — [zenodo 20768750](https://zenodo.org/records/20768750)
5. Hydrogeology (min–kyr) → super-GUE geyser, GOE glacial floods (**this work**)

## Honest Limitations
- Old Faithful dataset is 272 eruptions over ~2 weeks in August 1985; long-term drift not captured
- Bimodal waiting-time distribution complicates direct Wigner surmise comparison (⟨r⟩ remains valid)
- Multi-geyser model assumes independent Poisson components; real basins may have weak interactions
- IRD peak detection depends on threshold/prominence parameters (GOE classification robust across choices)
- IRD sample size modest (n=34); age model uncertainties ~0.5–2 kyr
- GOE-vs-GUE not resolved for IRD at present sample size

## Citation
```bibtex
@misc{chen2026hydro,
  author = {Chen, Ruqing},
  title  = {From Geysers to Megafloods: Random Matrix Theory Reveals
            Universal Level Repulsion in Hydrogeological
            Charge-and-Release Systems},
  year   = {2026},
  publisher = {GitHub},
  url    = {https://github.com/Ruqing1963/hydro-rmt-geysers-floods}
}
```

## License
[MIT](LICENSE). Old Faithful data from Azzalini & Bowman (1990).
DSDP 609 HSG data from PANGAEA (Obrochta et al. 2012, CC-BY-3.0).
