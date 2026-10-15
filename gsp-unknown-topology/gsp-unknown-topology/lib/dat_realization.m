function [fam, D] = dat_realization(cfg, caseName, famSeed, dataSeed, model, T, Ttest, famOpts)
%DAT_REALIZATION  One Monte Carlo realization: candidate family and data.
%
%   [fam, D] = DAT_REALIZATION(cfg, caseName, famSeed, dataSeed, model, T, Ttest)
%   draws the candidate family of caseName with seed famSeed (cfg.M topologies,
%   perturbation cfg.family unless famOpts is given) and, with seed dataSeed,
%   T training and Ttest test samples per topology from data model 'M1' or 'M2'.
%   Using the same famSeed in several experiments gives the same families.

if nargin < 8 || isempty(famOpts), famOpts = cfg.family; end
mpc = net_load_case(caseName);
ut_seed(famSeed);
fam = net_topology_family(mpc, cfg.M, famOpts);
ut_seed(dataSeed);
D = dat_family_data(fam, model, T, Ttest, cfg);
end
