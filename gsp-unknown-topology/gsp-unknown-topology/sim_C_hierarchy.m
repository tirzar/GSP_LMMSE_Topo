function res = sim_C_hierarchy(cfg)
%SIM_C_HIERARCHY  Section VI-C: number of bases and perturbation size (Figs. 2, 3).
%
%   res = SIM_C_HIERARCHY(cfg) runs two studies for data models M1 and M2 on
%   both networks:
%   (i)  number of bases: the multi-basis estimator (Theorem 1) restricted to
%        K of the M candidate bases (all singletons for K = 1, cfg.C.nSub
%        random subsets otherwise), compared with the topology-blind LMMSE
%        estimator (12), its limit by Theorem 2;
%   (ii) perturbation size: the candidate-family parameters (pRemove,
%        nAddMean) are scaled by cfg.C.scales, and the estimators cfg.C.keys
%        are evaluated. Scale 1 gives the families of Sections VI-B, VI-D, VI-E.
%   Output: results/sim_C.mat, figures/fig_C_bases, figures/fig_C_perturbation.

if nargin < 1, cfg = sim_config(); end
fprintf('\n=== Section VI-C: estimator hierarchy (%s mode) ===\n', cfg.mode);
models = {'M1', 'M2'};
snr = [cfg.snr0DB, cfg.snrM2DB];
nN = numel(cfg.networks); R = cfg.C.nMC; nK = numel(cfg.C.K); nSc = numel(cfg.C.scales);
nSub = max(cfg.C.nSub, cfg.M);
res = struct('nets', {cfg.networks}, 'netLabels', {cfg.netLabels}, 'models', {models}, ...
             'snrDB', snr, 'K', cfg.C.K, 'scales', cfg.C.scales, 'keys', {cfg.C.keys}, ...
             'M', cfg.M, 'T', cfg.T);
res.kb.nmse = NaN(nK, nSub, R, nN, 2); res.kb.worst = res.kb.nmse;
res.kb.ref = NaN(3, R, nN, 2); res.kb.refWorst = res.kb.ref;
res.kb.refKeys = {'blind_lmmse', 'expected_gsp', 'aware_gsp'};
nE = numel(cfg.C.keys);
res.pt.nmse = NaN(nE, nSc, R, nN, 2); res.pt.worst = res.pt.nmse;
res.pt.changes = NaN(nSc, R, nN); res.pt.align = res.pt.changes;
for mdl = 1:2
    sigma2 = ut_snr(snr(mdl));
    for n = 1:nN
        for r = 1:R
            t0 = tic;
            % (i) number of bases, default candidate family
            [fs, ds] = ut_seed_of(cfg, n, r, 3);
            [fam, D] = dat_realization(cfg, cfg.networks{n}, fs, ds, models{mdl}, cfg.T, cfg.Ttest);
            S = dat_moments(fam, D, sigma2, cfg.T);
            ests = est_design_set(res.kb.refKeys, fam, S, cfg);
            Rv = ev_family(ests, D, sigma2, fam.p);
            res.kb.ref(:, r, n, mdl) = Rv.nmse; res.kb.refWorst(:, r, n, mdl) = Rv.worst;
            ut_seed(ds + 7);
            for kk = 1:nK
                K = cfg.C.K(kk);
                if K == 1
                    subsets = num2cell(1:cfg.M);
                elseif K >= cfg.M
                    subsets = {1:cfg.M};
                else
                    subsets = cell(1, cfg.C.nSub);
                    for j = 1:cfg.C.nSub, subsets{j} = sort(randperm(cfg.M, K)); end
                end
                ests = cellfun(@(idx) est_multibasis(fam, S, idx), subsets, 'UniformOutput', false);
                Rv = ev_family(ests, D, sigma2, fam.p);
                res.kb.nmse(kk, 1:numel(subsets), r, n, mdl) = Rv.nmse;
                res.kb.worst(kk, 1:numel(subsets), r, n, mdl) = Rv.worst;
            end
            % (ii) perturbation size
            for sc = 1:nSc
                opts = cfg.family;
                opts.pRemove = opts.pRemove * cfg.C.scales(sc);
                opts.nAddMean = opts.nAddMean * cfg.C.scales(sc);
                variant = sc * (abs(cfg.C.scales(sc) - 1) > 1e-12);
                [fs, ds] = ut_seed_of(cfg, n, r, 3, variant);
                [fam, D] = dat_realization(cfg, cfg.networks{n}, fs, ds, models{mdl}, cfg.T, cfg.Ttest, opts);
                S = dat_moments(fam, D, sigma2, cfg.T);
                ut_seed(ds + 11);
                ests = est_design_set(cfg.C.keys, fam, S, cfg);
                Rv = ev_family(ests, D, sigma2, fam.p);
                res.pt.nmse(:, sc, r, n, mdl) = Rv.nmse; res.pt.worst(:, sc, r, n, mdl) = Rv.worst;
                res.pt.changes(sc, r, n) = mean(fam.nChanges(2:end));
                res.pt.align(sc, r, n) = mean(fam.align(2:end));
            end
            fprintf('  %s %-8s run %2d/%d (%.0f s)  multi-basis with K bases [dB]: %s| blind LMMSE %.2f\n', ...
                    models{mdl}, cfg.networks{n}, r, R, toc(t0), ...
                    sprintf('%.2f ', 10 * log10(mean(res.kb.nmse(:, :, r, n, mdl), 2, 'omitnan'))), ...
                    10 * log10(res.kb.ref(1, r, n, mdl)));
            save(fullfile(cfg.dirs.results, 'sim_C.mat'), 'res', '-v7');
        end
    end
end
printSummary(res);
fig_C_bases(res, cfg);
fig_C_perturbation(res, cfg);
end

function printSummary(res)
for mdl = 1:2
    for n = 1:numel(res.nets)
        fprintf('\n  %s, %s: multi-basis with K bases, mean / worst-topology NMSE [dB]\n', ...
                res.models{mdl}, res.netLabels{n});
        for kk = 1:numel(res.K)
            v = res.kb.nmse(kk, :, :, n, mdl); w = res.kb.worst(kk, :, :, n, mdl);
            fprintf('    K = %3d: %7.2f / %7.2f\n', res.K(kk), 10 * log10(mean(v(~isnan(v)))), ...
                    10 * log10(mean(w(~isnan(w)))));
        end
        for j = 1:numel(res.kb.refKeys)
            fprintf('    %-14s %7.2f / %7.2f\n', res.kb.refKeys{j}, ...
                    10 * log10(mean(res.kb.ref(j, :, n, mdl))), 10 * log10(mean(res.kb.refWorst(j, :, n, mdl))));
        end
        fprintf('  perturbation size, mean changed lines:%s\n', sprintf(' %6.2f', mean(res.pt.changes(:, :, n), 2)));
        for e = 1:numel(res.keys)
            fprintf('    %-16s%s\n', res.keys{e}, sprintf(' %6.2f', 10 * log10(mean(res.pt.nmse(e, :, :, n, mdl), 3))));
        end
    end
end
end
