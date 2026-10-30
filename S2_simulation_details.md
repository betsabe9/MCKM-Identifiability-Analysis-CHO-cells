# Supplementary S2 — Simulation details and comparison with the reference model

This note collects material moved from the main text for space reasons. Equation and figure numbers refer to the main text unless stated otherwise.

## S2.1 Modeling assumptions

The MCKM is a macroscopic, unstructured and unsegregated model. The following assumptions apply:

- The bioreactor is perfectly mixed, so temperature, pH and concentrations are spatially homogeneous.
- Cell growth is limited only by glucose, the primary carbon source. Glutamate is fed in excess and is not growth-limiting.
- Cell death is driven exclusively by the accumulation of cytotoxic ammonium in the medium.
- The culture is operated in fed-batch mode with variable volume, with discrete bolus feeding of glucose and glutamate over 15 days.

## S2.2 Numerical implementation

The bolus feeding couples continuous metabolic dynamics with discrete feed events, which makes the system stiff. The model was therefore simulated with CasADi (https://web.casadi.org/) and the implicit IDAS integrator, which uses backward differentiation formulas (BDF), with tolerances atol = 10⁻⁸ and rtol = 10⁻⁶.

**Table S1.** Initial conditions used in the simulations.

| X<sub>v</sub> | P | glc | glu | gln | amm | lac | μ |
|---|---|---|---|---|---|---|---|
| 0.3 × 10⁹ cells | 0.0 mM | 20.0 mM | 8.0 mM | 2.0 mM | 0.2 mM | 1.0 mM | 0.024 h⁻¹ |

All values were fixed except μ₀, which was computed from Eq. (2) with the initial glucose concentration.

## S2.3 Comparison with the reference simulations

**Fig. S2.** (a) Biomass and mAb profiles simulated in this work (CasADi-IDAS).

![Simulation in this work](FigS2a_this_work.png)

(b, c) Biomass and protein profiles reported by Sietaram et al. [1] (normalized units), reproduced from the ChemRxiv preprint under the CC BY 4.0 license.

![Sietaram et al. 2025, left panels](FigS2b_Sietaram_2025_left.png)

![Sietaram et al. 2025, right panels](FigS2c_Sietaram_2025_right.png)

The simulation reproduces the biomass growth and decline and the progressive mAb accumulation reported in [1].

## References

[1] D. Sietaram, P. Kotidis, G. Finka, and A. A. Lapkin, "Multi-cell line kinetic modelling of CHO fed-batch production runs to assess clone variability," ChemRxiv preprint, 2025. doi: 10.26434/chemrxiv-2025-pjhff. Published as "A multi clone kinetic model for characterizing Chinese hamster ovary cell line variability," *J. Ind. Microbiol. Biotechnol.*, vol. 52, 2025.
