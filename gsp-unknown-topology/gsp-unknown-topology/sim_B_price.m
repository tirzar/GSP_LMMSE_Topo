function res = sim_B_price(cfg)
%SIM_B_PRICE  Section VI-B: price of an unknown topology (Fig. 1, Table V).
%
%   res = SIM_B_PRICE(cfg) evaluates the estimators cfg.B.keys versus the SNR
%   on the IEEE 118-bus and New England 39-bus systems with M candidate
%   topologies (data model M1).
%   Output: results/sim_B.mat, figures/fig_B_price, tables/tab_B_comparison.tex.
%
%   res.nmse(e, s, r, n), res.worst, res.id, res.tDesign: estimator e
%   (res.keys), SNR s (res.snrDB), realization r, network n.

if nargin < 1, cfg = sim_config(); end
fprintf('\n=== Section VI-B: price of an unknown topology (%s mode) ===\n', cfg.mode);
keys = cfg.B.keys;
nE = numel(keys); nS = numel(cfg.snrDB); R = cfg.B.nMC; nN = numel(cfg.networks);
res = struct('keys', {keys}, 'snrDB', cfg.snrDB, 'nets', {cfg.networks}, ...
             'netLabels', {cfg.netLabels}, 'M', cfg.M, 'T', cfg.T, 'snr0DB', cfg.snr0DB);
res.nmse = NaN(nE, nS, R, nN); res.worst = res.nmse; res.id = res.nmse; res.tDesign = res.nmse;
for n = 1:nN
    for r = 1:R
        t0 = tic;
        [fs, ds] = ut_seed_of(cfg, n, r, 2);
        [fam, D] = dat_realization(cfg, cfg.networks{n}, fs, ds, 'M1', cfg.T, cfg.Ttest);
        for s = 1:nS
            sigma2 = ut_snr(cfg.snrDB(s));
            S = dat_moments(fam, D, sigma2, cfg.T);
            ut_seed(ds + s);
            [ests, td] = est_design_set(keys, fam, S, cfg);
            Rv = ev_family(ests, D, sigma2, fam.p);
            res.nmse(:, s, r, n) = Rv.nmse; res.worst(:, s, r, n) = Rv.worst;
            res.id(:, s, r, n) = Rv.id;     res.tDesign(:, s, r, n) = td;
        end
        s0 = cfg.snrDB == cfg.snr0DB;
        v = @(key) 10 * log10(res.nmse(strcmp(keys, key), s0, r, n));
        fprintf('  %-8s run %2d/%d (%.0f s)  NMSE at %g dB: aware %.2f, blind %.2f, MB %.2f, MM-GSP %.2f\n', ...
                cfg.networks{n}, r, R, toc(t0), cfg.snr0DB, v('aware_gsp'), v('blind_lmmse'), ...
                v('mb'), v('mm_gsp_soft'));
        save(fullfile(cfg.dirs.results, 'sim_B.mat'), 'res', '-v7');
    end
end
ut_print_summary(res, cfg.snr0DB);
fig_B_price(res, cfg);
tab_B_comparison(res, cfg);
end
