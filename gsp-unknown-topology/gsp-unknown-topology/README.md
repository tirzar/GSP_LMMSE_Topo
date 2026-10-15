# Graph Signal Estimation With Unknown Topology

MATLAB/GNU Octave code that reproduces the figures and tables of Section VI of

> T. Routtenberg, "Graph Signal Estimation With Unknown Topology," submitted to *IEEE Transactions on Signal Processing*, 2026.

Requires MATLAB R2019b or later (no toolboxes) or GNU Octave 8 or later. The IEEE test systems are included in `cases/`.

## Usage

```matlab
make_figures                        % figures and tables from the stored results (seconds)
run_all_simulations('quick')        % short test run with reduced sizes (minutes)
run_all_simulations                 % paper settings (several hours)
run_all_simulations('full', 'BE')   % selected parts only (here Sections VI-B and VI-E)
```

`results/` contains the results reported in the paper, so `make_figures` runs without new simulations. A full run overwrites them; it takes about 5 hours on one CPU core, and its parts can run in separate sessions. The quick mode uses fewer candidates, samples, and realizations, so its numbers differ from the paper. Figures are written to `figures/` (PDF and PNG) and LaTeX tables to `tables/`. All parameters are set in `sim_config.m`.

The stored results were computed with GNU Octave. MATLAB uses a different random number generator, so a MATLAB run gives statistically equivalent, but not identical, numbers.

## Contents

| Section | Script | Output | Paper |
|---|---|---|---|
| VI-A Setup | `sim_A_setup` | `tab_A_networks.tex` | Table III |
| VI-B Price of an unknown topology | `sim_B_price` | `fig_B_price`, `tab_B_comparison.tex` | Fig. 1, Table V |
| VI-C Number of bases and perturbation size | `sim_C_hierarchy` | `fig_C_bases`, `fig_C_perturbation` | Figs. 2, 3 |
| VI-D Finite training data | `sim_D_training` | `fig_D_training` | Fig. 4 |
| VI-E Topology information in the measurements | `sim_E_selection` | `fig_E_selection` | Fig. 5 |
| VI-F Transfer across network sizes | `sim_F_transfer` | `fig_F_collapse`, `tab_F_loo.tex` | Fig. 6, Table VI |

Each script prints a summary of its results, including the values quoted in the text (for example, the few-shot results of Section VI-F and the NE 39-bus results of Sections VI-D and VI-E), saves `results/sim_<X>.mat`, and calls its plotting and table functions (`fig_*`, `tab_*`).

## Estimators

| Key | Estimator | Reference in the paper |
|---|---|---|
| `mb` | Multi-basis GSP-LMMSE | Theorem 1 |
| `r1` | Rank-one | Remark 1 |
| `sr` | Shared response | Proposition 1 |
| `mm_gsp_hard` | Selective GSP (MAP topology selection) | (38) |
| `mm_gsp_soft` | Multiple-model GSP (posterior weighting) | (40) |
| `aware_gsp` | GSP-LMMSE, known topology (benchmark) | (33) |
| `aware_lmmse_ridge` | LMMSE, known topology (benchmark) | (14) |
| `blind_lmmse`, `blind_lmmse_ridge` | Topology-blind LMMSE | (12) |
| `nominal_gsp` | GSP-LMMSE of the nominal topology | Kroizer, Routtenberg, Eldar, 2022 |
| `expected_gsp` | Graph filter of the expected Laplacian | Isufi et al., 2017 |
| `minimax_mse`, `minimax_regret` | Minimax MSE and minimax regret | Eldar, Ben-Tal, Nemirovski, 2005; Eldar and Merhav, 2004 |
| `mm_full_hard`, `mm_full_soft` | Classical multiple-model estimator (full covariances) | Magill, 1965; Li and Jilkov, 2005 |
| `dnn` | Fully connected neural network | Mestav, Luengo-Rozas, Tong, 2019 |
| `sb_nom` | Single basis V<sub>1</sub> (ablation, Section VI-C) | Proposition 1 |

The implementation details of all estimators are given in the supplementary material of the paper.

## Code organization

| Folder / prefix | Content |
|---|---|
| `lib/net_*` | test cases, admittance matrix, Laplacian and graph Fourier basis, Laplacian checks, candidate topology families |
| `lib/dat_*` | data models M1 (smoothness prior) and M2 (AC power-flow states), sample moments |
| `lib/est_*` | design and application of all the estimators |
| `lib/ev_family.m` | NMSE, worst-topology NMSE, and topology-identification rate |
| `lib/xs_*` | cross-network design of Section IV: spectral dictionaries, fit, EVD-free filters |
| `lib/fg_*`, `lib/tab_write.m`, `lib/ut_*` | figure style, LaTeX tables, utilities |
| `cases/` | IEEE test systems in MATPOWER format |

## Notes

* SNR [dB] = −10 log<sub>10</sub> σ², with σ in p.u. on a 100-MVA base.
* Each Monte Carlo realization redraws the candidate family, the training data, and the test data. The seeds depend only on the network and the realization index, so Sections VI-B to VI-E use the same candidate families and a single realization can be re-run on its own.
* Curves show the Monte Carlo mean of the NMSE in dB. Error bars are 95% confidence intervals, ±1.96 standard errors of the per-realization NMSE in dB.
* Equation, section, figure, and table numbers in the code refer to the submitted manuscript.

## Test systems

The files in `cases/` are from MATPOWER (R. D. Zimmerman, C. E. Murillo-Sánchez, and R. J. Thomas, "MATPOWER: Steady-state operations, planning, and analysis tools for power systems research and education," *IEEE Trans. Power Syst.*, vol. 26, no. 1, pp. 12–19, 2011), distributed under the 3-clause BSD license.

## License

MIT, see `LICENSE`.
