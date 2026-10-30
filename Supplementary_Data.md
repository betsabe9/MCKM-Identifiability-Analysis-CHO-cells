# Supplementary S1 — Extended comparison between SVD-ROSM and the pairwise Collinearity Index

This note extends the discussion in Section III of the main text. It explains why the SVD of the time-resolved Relative Output Sensitivity Matrix (ROSM) [6] and the pairwise Collinearity Index (CI) used by Sietaram et al. [1] detect the same parameter pairs, how both relate to the general CI framework of Gábor et al. [4], and what additional information the SVD-ROSM framework provides. Equation numbers refer to the main text unless stated otherwise.

---

## S1.1 Relation between the CI and the SVD of the sensitivity matrix

For a subset `K` of parameters, the CI is defined as [2, 3]

<p align="center"><img src="https://latex.codecogs.com/png.image?%5Cdpi%7B200%7D%5Cbg%7Bwhite%7DCI_K%3D%5Cfrac%7B1%7D%7B%5Csqrt%7B%5Clambda_%7B%5Cmin%7D%5C%21%5Cleft%28%5Ctilde%7BS%7D_K%5E%7BT%7D%5Ctilde%7BS%7D_K%5Cright%29%7D%7D%3D%5Cfrac%7B1%7D%7B%5Csigma_%7B%5Cmin%7D%5C%21%5Cleft%28%5Ctilde%7BS%7D_K%5Cright%29%7D" alt="CI_K = 1/sqrt(lambda_min(S~_K^T S~_K)) = 1/sigma_min(S~_K)"></p>

where `S̃_K` is the sensitivity sub-matrix of the parameters in `K`, with columns normalized to unit length, and `σ_min` is its smallest singular value. The SVD of the full ROSM, `S = U Σ Vᵀ`, therefore evaluates the same quantity for the complete 13-parameter set, using relative rather than unit-norm column scaling. Both scalings remove the dependence on parameter units.

The CI is not restricted to pairs. Gábor et al. [4] evaluated it for parameter groups of any size, using combinatorial optimization to find the largest identifiable subset and incremental enumeration to find small groups of highly collinear parameters. In a TGF-β signalling model, they found a triplet of parameters with a very large CI even though none of its pairs was strongly collinear. This shows that pairwise screening alone can miss higher-order dependencies. In [1], the CI was applied only to parameter pairs.

Two consequences follow from the relation between the CI and the SVD:

- **Agreement is expected.** When a pair of parameters is collinear, both methods must detect it. The coincidence of the pairs identified here and in [1] therefore cross-validates both analyses rather than making one of them redundant.
- **The SVD also returns a direction.** The right singular vector `v_13` associated with `σ_13` gives the weights of the parameter combination that leaves the outputs unchanged. The CI only reports the degree of collinearity of a given subset.

## S1.2 Effect of time-averaging the sensitivities

In [1], the sensitivity of each output was averaged over the sampling times before computing the CI (Eq. 28 in [1]):

<p align="center"><img src="https://latex.codecogs.com/png.image?%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%5Cbar%7Bs%7D_%7Bij%7D%3D%5Cfrac%7B1%7D%7BT%7D%5Csum_%7Bt%3D1%7D%5E%7BT%7D%5Cfrac%7B%5Cpartial%20y_i%28t%29%7D%7B%5Cpartial%5Ctheta_j%7D" alt="s-bar_ij = (1/T) sum_t dy_i(t)/dtheta_j"></p>

which yields an `8 × 13` matrix (one row per output). In the present work, the ROSM keeps every sampling time as a separate row (`49 × 13`). Time-averaging is a choice of the implementation in [1], not a feature of the CI itself: in the formulation of Gábor et al. [4], each measured data point contributes a separate row of the sensitivity matrix, as in the ROSM.

In the MCKM, neither glutamine nor lactate appears in the balance of any other state: cell death depends only on ammonium (Eq. 4) and growth only on glucose (Eq. 2). Consequently:

- `Y_X/gln` and `Y_gln/glu` affect **only** the glutamine output;
- `Y_lac/glc` and `Y_X/lac` affect **only** the lactate output.

After time-averaging, each of these parameters has a sensitivity column with a single non-zero entry. Any two such columns that share the same non-zero row are proportional, so each pair is collinear **by construction**, independently of the underlying dynamics. This is consistent with the values reported in [1]: `CI = ∞` for the glutamine pair and `CI = 5×10⁶` for the lactate pair. The difference between both values is numerical, not structural.

Time-averaging also discards the information carried by the timing of each sensitivity. Sensitivities of opposite sign at different culture phases can partly cancel. More importantly, as shown below, the temporal signature is exactly what distinguishes a genuine structural degeneracy from an apparent one.

## S1.3 Glutamine: a genuine structural degeneracy

Substituting Eqs. (3) and (11) into Eq. (10) gives

<p align="center"><img src="https://latex.codecogs.com/png.image?%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%5Cfrac%7Bd%5C%2Cgln%7D%7Bdt%7D%3D%5Cmu%28t%29%5C%2CX_v%28t%29%5Cleft%28%5Cfrac%7BY_%7Bgln%2Fglu%7D%7D%7BY_%7BX%2Fglu%7D%7D-%5Cfrac%7B1%7D%7BY_%7BX%2Fgln%7D%7D%5Cright%29%5Cequiv%5Cmu%28t%29%5C%2CX_v%28t%29%5C%2CY_%7Bnet%7D" alt="d gln/dt = mu X_v (Y_gln/glu / Y_X/glu - 1/Y_X/gln) = mu X_v Y_net"></p>

The bracketed term is constant in time. `Y_X/glu` is determined independently by the glutamate balance (Eq. 8), so the glutamine output only carries information on the lumped quantity `Y_net`. Any change in `Y_gln/glu` can be exactly compensated by a change in `Y_X/gln` at every time point. This is a structural non-identifiability that no amount of noise-free data from the measured outputs can resolve.

**Biological interpretation.** `Y_gln/glu` describes glutamine synthesis from glutamate through glutamine synthetase (GS), whereas `Y_X/gln` describes glutamine consumption for biomass and protein synthesis. Both fluxes scale with growth (`μ·X_v`), so extracellular concentration data can only reveal their balance, i.e. the net glutamine accumulation per unit of growth. Separating synthesis from consumption would require additional information, such as isotopic tracers or intracellular flux measurements.

**Expected null-space direction.** In relative coordinates `ε_j = δθ_j / θ_j`, which are the coordinates of the ROSM, the condition `δY_net = 0` becomes

<p align="center"><img src="https://latex.codecogs.com/png.image?%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%5Cvarepsilon_%7Bgln%2Fglu%7D%5C%2C%5Cfrac%7BY_%7Bgln%2Fglu%7D%7D%7BY_%7BX%2Fglu%7D%7D%2B%5Cvarepsilon_%7BX%2Fgln%7D%5C%2C%5Cfrac%7B1%7D%7BY_%7BX%2Fgln%7D%7D%3D0" alt="eps_gln/glu * Y_gln/glu/Y_X/glu + eps_X/gln * 1/Y_X/gln = 0"></p>

With the nominal values of Table 1 (`Y_gln/glu = 0.38`, `Y_X/glu = 8.9`, `Y_X/gln = 43.2`):

<p align="center"><img src="https://latex.codecogs.com/png.image?%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%5Cfrac%7BY_%7Bgln%2Fglu%7D%7D%7BY_%7BX%2Fglu%7D%7D%5Capprox0.0427%2C%5Cqquad%5Cfrac%7B1%7D%7BY_%7BX%2Fgln%7D%7D%5Capprox0.0231%2C%5Cqquad%5Cfrac%7B%5Cvarepsilon_%7BX%2Fgln%7D%7D%7B%5Cvarepsilon_%7Bgln%2Fglu%7D%7D%5Capprox-1.85" alt="Y_gln/glu / Y_X/glu ~ 0.0427; 1/Y_X/gln ~ 0.0231; eps_X/gln / eps_gln/glu ~ -1.85"></p>

The normalized null-space vector is therefore expected to be approximately

<p align="center"><img src="https://latex.codecogs.com/png.image?%5Cdpi%7B200%7D%5Cbg%7Bwhite%7Dv_%7B13%7D%5Capprox%5Cpm%5Cleft%280.48%5C%20%5Ctext%7Bfor%7D%5C%20Y_%7Bgln%2Fglu%7D%2C%5C%20-0.88%5C%20%5Ctext%7Bfor%7D%5C%20Y_%7BX%2Fgln%7D%5Cright%29" alt="v_13 ~ +/-(0.48 for Y_gln/glu, -0.88 for Y_X/gln)"></p>

with zero weight on the remaining parameters. The opposite signs mean that a higher glutamine synthesis is compensated by a higher glutamine consumption. The numerical null-space vector shown in Fig. 2 of the main text can be compared directly with this analytical prediction.

Note also that `Y_net ≈ 0.0427 − 0.0231 = 0.0196 > 0`, i.e. net glutamine accumulation. This agrees with the increasing glutamine profiles reported for GS-CHO cultures in [1].

**Model reduction.** Fixing `Y_gln/glu`, whose value is bounded above by the 1:1 stoichiometry of the GS reaction, removes this direction and leaves `Y_X/gln` as the clone-specific parameter. An equivalent alternative would be to reparameterize the model and estimate `Y_net` directly.

## S1.4 Lactate: structurally identifiable, practically limited

From Eqs. (14)–(16), the specific lactate rate has two regimes:

<p align="center"><img src="https://latex.codecogs.com/png.image?%5Cdpi%7B200%7D%5Cbg%7Bwhite%7Dq_%7Blac%7D%28t%29%3D%5Cbegin%7Bcases%7DY_%7Blac%2Fglc%7D%5C%2Cq_%7Bglc%7D%28t%29%2C%20%26%20glc%28t%29%3E%5Calpha_%7Bglc%7D%5Cquad%5Ctext%7B%28production%29%7D%5C%5C%5B6pt%5D%20-%5Cdfrac%7B%5Cmu%28t%29%7D%7BY_%7BX%2Flac%7D%7D%2BY_%7Blac%2Fglc%7D%5C%2Cq_%7Bglc%7D%28t%29%2C%20%26%20glc%28t%29%5Cle%5Calpha_%7Bglc%7D%5Cquad%5Ctext%7B%28after%20the%20switch%29%7D%5Cend%7Bcases%7D" alt="q_lac(t), two regimes (production / after the switch)"></p>

The sensitivity of lactate to `Y_lac/glc` is non-zero throughout the culture. The sensitivity to `Y_X/lac` is zero before the switch and only builds up afterwards. Their columns in the time-resolved ROSM are therefore **not** proportional. In other words, `Y_lac/glc` is informed by the production phase, and once it is known, `Y_X/lac` is informed by the consumption phase. This is why the SVD does not flag this pair as structurally non-identifiable, even though it appears strongly collinear in the time-averaged analysis of [1].

Two qualifications apply:

1. **Run-dependent identifiability.** Structural identifiability of `Y_X/lac` requires the switch to occur within the sampled window. For clones that only produce lactate (as shown in [1], Fig. 3), the sensitivity column of `Y_X/lac` is zero and the parameter is structurally non-identifiable for that run.
2. **Practical limitation.** After the switch, glucose is low, so `μ` and hence the term `μ / Y_X/lac` are small. In addition, only the few samples taken after glucose depletion (out of seven, at days 0, 3, 6, 8, 10, 13 and 15) carry information on lactate consumption. This explains the low t-value of `Y_X/lac` under 5 % measurement noise in the main text.

## S1.5 Hierarchy and computational cost

**One degeneracy, one fixed parameter.** For the full model, only one singular value collapses to machine precision (`σ_13 ≈ 1×10⁻¹⁵`). Removing one column from a matrix interlaces its singular values [5]:

<p align="center"><img src="https://latex.codecogs.com/png.image?%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%5Csigma_k%28S%29%5C%20%5Cge%5C%20%5Csigma_k%28S_%7Bred%7D%29%5C%20%5Cge%5C%20%5Csigma_%7Bk%2B1%7D%28S%29" alt="sigma_k(S) >= sigma_k(S_red) >= sigma_(k+1)(S)"></p>

Therefore `σ_12(S) ≥ σ_12(S_red) ≈ 5×10⁻³`, so the spectral gap between `σ_12` and `σ_13` spans at least twelve orders of magnitude. The full model has exactly one structural degeneracy, and fixing a single parameter in the dominant direction removes it. Because the gap is so large, this conclusion does not depend on the particular threshold chosen (`10⁻⁴` in the main text).

**Pairwise CI.** The pairwise CI yields `C(13,2) = 78` independent values. It does not indicate whether fixing one parameter resolves other dependencies, nor which parameter in a pair should be fixed. Higher-order dependencies can be detected with the CI by evaluating larger subsets [4]. For the 13 parameters of the MCKM, all `2¹³ − 13 − 1 = 8178` subsets of two or more parameters could still be enumerated exhaustively.

The advantage of SVD-ROSM in this model is therefore not computational speed. It is that a single decomposition returns the non-identifiable direction together with the weight of each parameter (S1.3), which a CI value does not provide. For larger models, exhaustive enumeration of collinear groups becomes infeasible beyond a few dozen parameters [4], whereas the cost of an SVD grows only polynomially with the number of parameters.

## S1.6 Mechanistic interpretation of practical non-identifiability

The t-value analysis of the main text identified `K_glc`, `Y_X/lac` and `m_glc` as practically non-identifiable in the 12-parameter model. `Y_X/lac` is discussed in S1.4. The other two have direct mechanistic explanations.

**`K_glc`.** The relative sensitivity of the growth rate (Eq. 2) to `K_glc` is

<p align="center"><img src="https://latex.codecogs.com/png.image?%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%5Cfrac%7BK_%7Bglc%7D%7D%7B%5Cmu%7D%5Cfrac%7B%5Cpartial%5Cmu%7D%7B%5Cpartial%20K_%7Bglc%7D%7D%3D-%5Cfrac%7BK_%7Bglc%7D%7D%7BK_%7Bglc%7D%2Bglc%7D" alt="(K_glc/mu) dmu/dK_glc = -K_glc/(K_glc+glc)"></p>

With `glc ≈ 20 mM` and `K_glc ≈ 1.05 mM`, this is about `−0.05`. It only approaches `−0.5` when glucose falls to the order of `K_glc`, which under bolus feeding happens during short periods between feeds. While glucose is high, the small and nearly constant effect of `K_glc` can be absorbed by `μ_max`, whose relative sensitivity on `μ` is always 1. This explains both the weak practical identifiability of `K_glc` and the `μ_max`–`K_glc` correlation reported in [1] (`CI = 14.33`).

**`m_glc`.** In Eq. (9), the maintenance term represents a fraction

<p align="center"><img src="https://latex.codecogs.com/png.image?%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%5Cfrac%7Bm_%7Bglc%7D%7D%7B%5Cmu%2FY_%7BX%2Fglc%7D%2Bm_%7Bglc%7D%7D%5Capprox%5Cfrac%7B3.15%5Ctimes10%5E%7B-6%7D%7D%7B0.0127%7D%5Capprox2.5%5Ctimes10%5E%7B-4%7D" alt="m_glc/(mu/Y_X/glc + m_glc) ~ 3.15e-6/0.0127 ~ 2.5e-4"></p>

of the specific glucose uptake during growth (nominal values of Table 1). Its influence on the outputs only becomes appreciable when `μ` approaches zero, so the data contain very little information on this parameter.

**Implications for experimental design.** The three practically non-identifiable parameters share the same information window: the periods around glucose depletion and after the lactate switch. Denser sampling in that window (for example, between days 6 and 10) is the most direct route to make them estimable without fixing them. By contrast, the glutamine degeneracy (S1.3) cannot be resolved by more frequent sampling of the same outputs.

Poor practical identifiability is common in kinetic models of CHO metabolism. In a large-scale CHO benchmark model with more than 100 parameters, Gábor et al. [4] found that only 13 parameters were practically identifiable with the available data.

## S1.7 Implications for inter-clonal comparison

The purpose of the MCKM is to compare clones through their kinetic parameters. In [1], linear discriminant analysis (LDA) of the regressed parameters was used to separate high- and low-performing clones and stable and unstable clones (Fig. 5 in [1]). Several parameters with substantial LDA weight, such as `K_glc`, are practically non-identifiable according to the present analysis, and `Y_X/gln` is only identifiable once `Y_gln/glu` is fixed. Differences between clones in these parameters may therefore partly reflect estimation uncertainty rather than metabolic differences.

We suggest restricting inter-clonal comparisons to the identifiable subset of the 9-parameter model (`μ_max`, `k_d`, `Y_X/glc`, `Y_X/glu`, `Y_X/gln`, `Y_lac/glc`, `Y_X/amm`, `Y_P/X`, `KD_amm`), or to identifiable combinations such as `Y_net`.

## S1.8 Limitations

- The SVD-ROSM analysis is local: it is evaluated around nominal parameter values. The Monte Carlo analysis (±20 %, 100 sets) supports the robustness of the structural conclusions within that range, but not globally.
- Structural conclusions assume that all seven outputs are measured at the seven sampling times of [1].
- The practical identifiability analysis assumes a uniform relative measurement error of 5 %. Different noise levels, or output-specific errors, may change which parameters fall below the t-value threshold.
- As noted in S1.4, the identifiability of `Y_X/lac` depends on whether a given clone exhibits the lactate switch.

---

## References

[1] D. Sietaram, P. Kotidis, G. Finka, and A. A. Lapkin, "Multi-cell line kinetic modelling of CHO fed-batch production runs to assess clone variability," ChemRxiv preprint, 2025. doi: 10.26434/chemrxiv-2025-pjhff.

[2] R. Brun, P. Reichert, and H. R. Künsch, "Practical identifiability analysis of large environmental simulation models," *Water Resources Research*, vol. 37, no. 4, pp. 1015–1030, 2001.

[3] R. Brun, M. Kühni, H. Siegrist, W. Gujer, and P. Reichert, "Practical identifiability of ASM2d parameters—systematic selection and tuning of parameter subsets," *Water Research*, vol. 36, no. 16, pp. 4113–4127, 2002. doi: 10.1016/S0043-1354(02)00104-5.

[4] A. Gábor, A. F. Villaverde, and J. R. Banga, "Parameter identifiability analysis and visualization in large-scale kinetic models of biosystems," *BMC Systems Biology*, vol. 11, art. 54, 2017. doi: 10.1186/s12918-017-0428-y.

[5] R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed. Cambridge, U.K.: Cambridge University Press, 2013.

[6] J. D. Stigter and J. Molenaar, "A fast algorithm to assess local structural identifiability," *Automatica*, vol. 58, pp. 118–124, 2015.
