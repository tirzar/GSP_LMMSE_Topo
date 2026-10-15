function res = sim_E_selection(cfg)
%SIM_E_SELECTION  Section VI-E: topology information in the measurements (Fig. 5).
%
%   res = SIM_E_SELECTION(cfg) evaluates the selective (hard) and
%   multiple-model (soft) estimators when P snapshots share the active
%   topology, for P in cfg.E.P and T_m in cfg.E.T training samples per
%   topology (data model M1, SNR cfg.snr0DB). The log-likelihoods of the P
%   snapshots are summed. The graph-stationary surrogate likelihood (GSP-LMMSE
%   branches) is compared with the classical multiple-model estimator (full
%   covariance likelihoods and LMMSE branches).
%   Output: results/sim_E.mat, figures/fig_E_selection.
%
%   res.nmse(e, k, t, r, n), res.id(...): estimator e, number of snapshots
%   res.P(k), training size res.T(t), realization r, network n.

if nargin < 1, cfg = sim_config(); end
fprintf('\n=== Section VI-E: topology selection from P snapshots (%s mode) ===\n', cfg.mode);
keys = cfg.E.keys; nE = numel(keys); Ps = cfg.E.P; nP = numel(Ps);
Ts = cfg.E.T; nT = numel(Ts); nN = numel(cfg.networks); R = cfg.E.nMC;
sigma2 = ut_snr(cfg.snr0DB);
res = struct('nets', {cfg.networks}, 'netLabels', {cfg.netLabels}, 'keys', {keys}, ...
             'P', Ps, 'T', Ts, 'snrDB', cfg.snr0DB, 'M', cfg.M);
res.nmse = NaN(nE, nP, nT, R, nN); res.worst = res.nmse; res.id = res.nmse;
for n = 1:nN
    for r = 1:R
        t0 = tic;
        [fs, ds] = ut_seed_of(cfg, n, r, 5);
        [fam, D] = dat_realization(cfg, cfg.networks{n}, fs, ds, 'M1', max(Ts), cfg.E.Ttest);
        for t = 1:nT
            S = dat_moments(fam, D, sigma2, Ts(t));
            ut_seed(ds + t);
            ests = est_design_set(keys, fam, S, cfg);
            for k = 1:nP
                Rv = ev_family(ests, D, sigma2, fam.p, Ps(k));
                res.nmse(:, k, t, r, n) = Rv.nmse; res.worst(:, k, t, r, n) = Rv.worst;
                res.id(:, k, t, r, n) = Rv.id;
            end
        end
        e = strcmp(keys, 'mm_gsp_hard');
        fprintf('  %-8s run %2d/%d (%.0f s)  identification rate (T_m = %d): %s\n', cfg.networks{n}, ...
                r, R, toc(t0), Ts(end), sprintf('%.2f ', res.id(e, :, end, r, n)));
        save(fullfile(cfg.dirs.results, 'sim_E.mat'), 'res', '-v7');
    end
end
for n = 1:nN
    for t = 1:nT
        fprintf('\n  %s, T_m = %d: NMSE [dB] | identification rate, P = %s\n', cfg.netLabels{n}, ...
                Ts(t), mat2str(Ps));
        for e = 1:nE
            V = reshape(res.nmse(e, :, t, :, n), nP, []);
            I = mean(reshape(res.id(e, :, t, :, n), nP, []), 2)';
            fprintf('    %-14s%s |%s\n', keys{e}, sprintf(' %7.2f', ut_mc_rows(V)), sprintf(' %5.2f', I));
        end
    end
end
fig_E_selection(res, cfg);
end
