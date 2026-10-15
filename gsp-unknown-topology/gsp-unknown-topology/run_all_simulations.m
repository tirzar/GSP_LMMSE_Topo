function run_all_simulations(mode, parts)
%RUN_ALL_SIMULATIONS  Reproduce the numerical results of Section VI.
%
%   RUN_ALL_SIMULATIONS()              paper settings (several hours)
%   RUN_ALL_SIMULATIONS('quick')       reduced sizes (a few minutes)
%   RUN_ALL_SIMULATIONS(mode, 'BD')    only Sections VI-B and VI-D
%
%   Each part saves results/sim_<X>.mat and writes its figures and tables:
%     A  SIM_A_SETUP       test systems and candidate families     Table III
%     B  SIM_B_PRICE       NMSE versus SNR, estimator comparison   Fig. 1, Table V
%     C  SIM_C_HIERARCHY   number of bases, perturbation size      Figs. 2, 3
%     D  SIM_D_TRAINING    finite training data                    Fig. 4
%     E  SIM_E_SELECTION   P snapshots sharing the topology        Fig. 5
%     F  SIM_F_TRANSFER    transfer across network sizes           Fig. 6, Table VI

if nargin < 1 || isempty(mode), mode = 'full'; end
if nargin < 2 || isempty(parts), parts = 'ABCDEF'; end
cfg = sim_config(mode);
fun = struct('A', @sim_A_setup, 'B', @sim_B_price, 'C', @sim_C_hierarchy, ...
             'D', @sim_D_training, 'E', @sim_E_selection, 'F', @sim_F_transfer);
tAll = tic;
for p = upper(parts)
    t0 = tic;
    feval(fun.(p), cfg);
    fprintf('--- Section VI-%s done in %.1f min\n', p, toc(t0) / 60);
end
fprintf('=== done in %.1f min; figures in %s, tables in %s\n', toc(tAll) / 60, ...
        cfg.dirs.figures, cfg.dirs.tables);
end
