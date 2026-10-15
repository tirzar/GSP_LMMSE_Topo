function make_figures(mode)
%MAKE_FIGURES  Re-create the figures and tables from the saved results.
%
%   MAKE_FIGURES() plots from results/sim_*.mat without re-running the
%   simulations. Missing result files are skipped. Figures are written to
%   figures/*.pdf and *.png, tables to tables/*.tex.

if nargin < 1, mode = 'full'; end
cfg = sim_config(mode);
f = @(x) fullfile(cfg.dirs.results, x);
if exist(f('sim_A.mat'), 'file'), S = load(f('sim_A.mat')); tab_A_networks(S.res, cfg); end
if exist(f('sim_B.mat'), 'file')
    S = load(f('sim_B.mat')); fig_B_price(S.res, cfg); tab_B_comparison(S.res, cfg);
end
if exist(f('sim_C.mat'), 'file')
    S = load(f('sim_C.mat')); fig_C_bases(S.res, cfg); fig_C_perturbation(S.res, cfg);
end
if exist(f('sim_D.mat'), 'file'), S = load(f('sim_D.mat')); fig_D_training(S.res, cfg); end
if exist(f('sim_E.mat'), 'file'), S = load(f('sim_E.mat')); fig_E_selection(S.res, cfg); end
if exist(f('sim_F.mat'), 'file')
    S = load(f('sim_F.mat')); fig_F_collapse(S.res, cfg); tab_F_loo(S.res, cfg);
end
end
