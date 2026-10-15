function res = sim_D_training(cfg)
%SIM_D_TRAINING  Section VI-D: finite training data (Fig. 4).
%
%   res = SIM_D_TRAINING(cfg) evaluates the estimators cfg.D.keys versus the
%   number of training samples per topology T_m (cfg.D.T), for data models M1
%   (SNR cfg.snr0DB) and M2 (SNR cfg.snrM2DB) on both networks. The training
%   sets are nested (the first T_m samples of one draw). The DNN is trained
%   only for T_m in cfg.D.dnnT.
%   Output: results/sim_D.mat, figures/fig_D_training.
%
%   res.nmse(e, t, r, n, mdl): estimator e, training size t, realization r,
%   network n, data model mdl.

if nargin < 1, cfg = sim_config(); end
fprintf('\n=== Section VI-D: finite training data (%s mode) ===\n', cfg.mode);
models = {'M1', 'M2'};
snr = [cfg.snr0DB, cfg.snrM2DB];
keys = cfg.D.keys; nE = numel(keys); Ts = cfg.D.T; nT = numel(Ts);
nN = numel(cfg.networks); R = cfg.D.nMC;
res = struct('nets', {cfg.networks}, 'netLabels', {cfg.netLabels}, 'models', {models}, ...
             'snrDB', snr, 'T', Ts, 'keys', {keys}, 'M', cfg.M);
res.nmse = NaN(nE, nT, R, nN, 2); res.worst = res.nmse; res.tDesign = res.nmse;
iDnn = strcmp(keys, 'dnn');
for mdl = 1:2
    sigma2 = ut_snr(snr(mdl));
    for n = 1:nN
        for r = 1:R
            t0 = tic;
            [fs, ds] = ut_seed_of(cfg, n, r, 4);
            [fam, D] = dat_realization(cfg, cfg.networks{n}, fs, ds, models{mdl}, max(Ts), cfg.Ttest);
            for t = 1:nT
                S = dat_moments(fam, D, sigma2, Ts(t));
                k = ~iDnn | any(Ts(t) == cfg.D.dnnT);
                ut_seed(ds + t);
                [ests, td] = est_design_set(keys(k), fam, S, cfg);
                Rv = ev_family(ests, D, sigma2, fam.p);
                res.nmse(k, t, r, n, mdl) = Rv.nmse; res.worst(k, t, r, n, mdl) = Rv.worst;
                res.tDesign(k, t, r, n, mdl) = td;
            end
            fprintf('  %s %-8s run %2d/%d (%.0f s)  multi-basis [dB]: %s\n', models{mdl}, ...
                    cfg.networks{n}, r, R, toc(t0), ...
                    sprintf('%.2f ', 10 * log10(res.nmse(strcmp(keys, 'mb'), :, r, n, mdl))));
            save(fullfile(cfg.dirs.results, 'sim_D.mat'), 'res', '-v7');
        end
    end
end
for mdl = 1:2
    for n = 1:nN
        fprintf('\n  %s, %s: NMSE [dB] versus T_m = %s\n', models{mdl}, cfg.netLabels{n}, mat2str(Ts));
        for e = 1:nE
            V = reshape(res.nmse(e, :, :, n, mdl), nT, []);
            fprintf('    %-18s%s\n', keys{e}, sprintf(' %7.2f', ut_mc_rows(V)));
        end
    end
end
fig_D_training(res, cfg);
end
