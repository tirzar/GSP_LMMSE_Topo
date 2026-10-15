# Graph Signal Estimation With Unknown Topology

MATLAB/GNU Octave code for reproducing the figures and tables in Section VI of:

> T. Routtenberg, "Graph Signal Estimation With Unknown Topology," submitted to *IEEE Transactions on Signal Processing*, 2026.

The code implements graph-signal estimators under topology uncertainty, including multi-basis, shared-response, low-rank, and topology-adaptive estimators, as well as shared spectral-function designs for transfer across networks of different sizes.

**Requirements:** MATLAB R2019b or later (no additional toolboxes) or GNU Octave 8 or later. The IEEE test systems are included in `cases/`.

## Usage

```matlab
make_figures                        % Reproduce figures and tables from stored results (seconds)
run_all_simulations('quick')        % Short test run with reduced settings (minutes)
run_all_simulations                 % Full paper settings (several hours)
run_all_simulations('full', 'BE')   % Run selected parts (Sections VI-B and VI-E)
```

The `results/` directory contains the numerical results reported in the paper, allowing `make_figures` to reproduce the figures and tables without rerunning the simulations.

A full simulation run overwrites the stored results and takes approximately 5 hours on a single CPU core. Individual parts can be executed in separate sessions. Quick mode uses fewer candidate topologies, samples, and Monte Carlo realizations; therefore, its numerical results differ from those reported in the paper.

Figures are saved to `figures/` in PDF and PNG formats, and LaTeX tables are saved to `tables/`. Simulation parameters are specified in `sim_config.m`.

The stored results were generated using GNU Octave. MATLAB uses a different random number generator, so MATLAB runs are expected to produce statistically comparable, but not identical, results.

## Reproducing the Paper

| Section | Simulation script | Output | Paper |
|---|---|---|---|
| VI-A: Setup | `sim_A_setup` | `tab_A_networks.tex` | Table III |
| VI-B: Price of an unknown topology | `sim_B_price` | `fig_B_price`, `tab_B_comparison.tex` | Fig. 1, Table V |
| VI-C: Number of bases and perturbation size | `sim_C_hierarchy` | `fig_C_bases`, `fig_C_perturbation` | Figs. 2–3 |
| VI-D: Finite training data | `sim_D_training` | `fig_D_training` | Fig. 4 |
| VI-E: Topology information in measurements | `sim_E_selection` | `fig_E_selection` | Fig. 5 |
| VI-F: Cross-network transfer | `sim_F_transfer` | `fig_F_collapse`, `tab_F_loo.tex` | Fig. 6, Table VI |

Each simulation script prints a summary of its results, saves a corresponding `
