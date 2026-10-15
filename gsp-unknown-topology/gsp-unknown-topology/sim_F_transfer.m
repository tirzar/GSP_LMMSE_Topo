function res = sim_F_transfer(cfg)
%SIM_F_TRANSFER  Section VI-F: transfer across network sizes (Fig. 6, Table VI).
%
%   res = SIM_F_TRANSFER(cfg) uses the six networks cfg.F.networks, each with
%   its known topology (one graph per group, Remark 4), data models M1 (SNR
%   cfg.F.snrDB) and M2 (SNR cfg.snrM2DB), and cfg.F.nMC independent data sets:
%   (1) GSP-LMMSE responses (55) of all networks versus the graph frequency;
%   (2) leave-one-network-out transfer: a spectral function is fitted on five
%       networks (Theorem 3, (56)) and applied to the sixth by (57), with
%         'cheb_abs'      Chebyshev, absolute frequency axis
%         'cheb_lmax'     Chebyshev, axis lambda / lambda_max (ChebNet)
%         'cheb_graphon'  Chebyshev, axis lambda / N (graphon)
%         'rational'      rational dictionary with fixed poles
%         'wiener'        beta / (beta lambda + sigma^2), linearized M1 (no learning)
%       The excess MSE over the target's own GSP-LMMSE estimator is measured on
%       test data with the EVD-free filters, and predicted from the training
%       moments by completing the square;
%   (3) few-shot adaptation with T0 labeled samples of the held-out network:
%         'target_gsp'   GSP-LMMSE estimator fitted on the T0 samples
%         'target_rat'   rational function fitted on the T0 samples
%         'transfer'     rational function learned on the other five networks
%         'joint'        rational function fitted jointly on the five networks
%                        and the T0 target samples (Section V-B)
%   Output: results/sim_F.mat, figures/fig_F_collapse, tables/tab_F_loo.tex.

if nargin < 1, cfg = sim_config(); end
fprintf('\n=== Section VI-F: transfer across network sizes (%s mode) ===\n', cfg.mode);
F = cfg.F;
J = numel(F.networks); R = F.nMC;
models = {'M1', 'M2'};
snr = [F.snrDB, cfg.snrM2DB];
dictNames = {'cheb_abs', 'cheb_lmax', 'cheb_graphon', 'rational', 'wiener'};
fewNames = {'target_gsp', 'target_rat', 'transfer', 'joint'};
nD = numel(dictNames); nF = numel(fewNames); nT0 = numel(F.T0);
res = struct('nets', {F.networks}, 'netLabels', {F.netLabels}, 'models', {models}, 'snrDB', snr, ...
             'dictNames', {dictNames}, 'fewNames', {fewNames}, 'T0', F.T0, 'Q', F.Q, 'beta', cfg.beta);
res.N = zeros(1, J); res.lamMax = zeros(1, J);
res.collapse = cell(J, 2);                     % [lambda, f_1, ..., f_R] per network and model
res.loo.own = NaN(J, R, 2);                    % test MSE of the target's own GSP-LMMSE estimator
res.loo.den = NaN(J, R, 2);                    % E||x - mu_x||^2 (test)
res.loo.meas = NaN(J, nD, R, 2);               % test MSE, EVD-free filters
res.loo.pred = NaN(J, nD, R, 2);               % predicted from the training moments
res.loo.evd = NaN(J, nD, R, 2);                % test MSE, with the EVD
res.few.mse = NaN(J, nF, nT0, F.nFew, R, 2);
for r = 1:R
    for mdl = 1:2
        t0 = tic;
        sigma2 = ut_snr(snr(mdl));
        nets = cell(1, J); tr = cell(1, J); te = cell(1, J); lam = cell(1, J);
        for j = 1:J
            nets{j} = xs_network_data(cfg, F.networks{j}, models{mdl}, F.T, F.Ttest, ...
                                      cfg.seed + 7000000 + 100 * r + 10 * mdl + j);
            g = nets{j}.grd;
            tr{j} = xs_moments(g, nets{j}.X, nets{j}.G + sqrt(sigma2) * nets{j}.W);
            te{j} = xs_moments(g, nets{j}.Xte, nets{j}.Gte + sqrt(sigma2) * nets{j}.Wte);
            lam{j} = g.lambda;
            res.N(j) = g.N; res.lamMax(j) = g.lambda(end);
            if r == 1
                res.collapse{j, mdl} = [g.lambda, tr{j}.f];
            else
                res.collapse{j, mdl} = [res.collapse{j, mdl}, tr{j}.f];
            end
            res.loo.den(j, r, mdl) = te{j}.trCxx;
            res.loo.own(j, r, mdl) = xs_eval(te{j}, g.V, tr{j}.mux, tr{j}.muy, tr{j}.f);
        end
        dicts = makeDicts(F, res.lamMax, res.N, sigma2, cfg.beta);
        % (2) leave-one-network-out
        alphaXS = cell(J, nD);
        for j0 = 1:J
            train = setdiff(1:J, j0);
            g = nets{j0}.grd;
            Yte = nets{j0}.Gte + sqrt(sigma2) * nets{j0}.Wte;
            for d = 1:nD
                if strcmp(dictNames{d}, 'wiener')
                    if mdl == 2, continue; end
                    alpha = 1;
                else
                    alpha = xs_fit(tr, lam, train, dicts{d});
                end
                alphaXS{j0, d} = alpha;
                h = xs_response(dicts{d}, alpha, lam{j0}, j0);
                res.loo.pred(j0, d, r, mdl) = tr{j0}.mseGSP + ...
                    sum(tr{j0}.dyy(2:end) .* (h(2:end) - tr{j0}.f(2:end)) .^ 2);
                [res.loo.evd(j0, d, r, mdl), res.loo.meas(j0, d, r, mdl)] = ...
                    xs_test_mse(g, tr{j0}, h, nets{j0}.Xte, Yte, dicts{d}, alpha, j0);
            end
        end
        % (3) few-shot adaptation, rational dictionary
        iRat = strcmp(dictNames, 'rational');
        dR = dicts{iRat};
        ut_seed(cfg.seed + 8000000 + 100 * r + mdl);
        for j0 = 1:J
            train = setdiff(1:J, j0);
            g = nets{j0}.grd;
            for k = 1:nT0
                for b = 1:F.nFew
                    idx = randperm(F.T, F.T0(k));
                    t0r = xs_moments(g, nets{j0}.X(:, idx), nets{j0}.G(:, idx) + sqrt(sigma2) * nets{j0}.W(:, idx));
                    res.few.mse(j0, 1, k, b, r, mdl) = xs_eval(te{j0}, g.V, t0r.mux, t0r.muy, t0r.f);
                    a = xs_fit({t0r}, {lam{j0}}, 1, dR);
                    res.few.mse(j0, 2, k, b, r, mdl) = xs_eval(te{j0}, g.V, t0r.mux, t0r.muy, ...
                                                               xs_response(dR, a, lam{j0}, 1));
                    h = xs_response(dR, alphaXS{j0, iRat}, lam{j0}, 1);
                    res.few.mse(j0, 3, k, b, r, mdl) = xs_eval(te{j0}, g.V, t0r.mux, t0r.muy, h);
                    a = xs_fit([tr(train), {t0r}], [lam(train), {lam{j0}}], 1:J, dR);
                    res.few.mse(j0, 4, k, b, r, mdl) = xs_eval(te{j0}, g.V, t0r.mux, t0r.muy, ...
                                                               xs_response(dR, a, lam{j0}, 1));
                end
            end
        end
        ex = 10 * log10(res.loo.meas(:, :, r, mdl) ./ res.loo.own(:, r, mdl));
        fprintf('  %s data set %d/%d (%.0f s)  excess MSE [dB] for %s:\n', models{mdl}, r, R, ...
                toc(t0), strjoin(F.networks, ', '));
        for d = 1:nD
            fprintf('    %-13s%s\n', dictNames{d}, sprintf(' %7.2f', ex(:, d)));
        end
        save(fullfile(cfg.dirs.results, 'sim_F.mat'), 'res', '-v7');
    end
end
printSummary(res);
fig_F_collapse(res, cfg);
tab_F_loo(res, cfg);
end

function dicts = makeDicts(F, lamMax, N, sigma2, beta)
% Dictionaries in the order of dictNames.
sAbs = max(lamMax) / 2;                     % common absolute axis covering all spectra
sLmax = lamMax / 2;                         % lambda / lambda_max
sGraphon = N * max(lamMax ./ N) / 2;        % lambda / N
dicts = {xs_dictionary('chebyshev', F.Q, sAbs, F.muCheb), ...
         xs_dictionary('chebyshev', F.Q, sLmax, F.muCheb), ...
         xs_dictionary('chebyshev', F.Q, sGraphon, F.muCheb), ...
         xs_dictionary('rational', F.poles, F.muRat), ...
         xs_dictionary('rational', sigma2 / beta, 0)};
end

function printSummary(res)
for mdl = 1:2
    fprintf('\n  %s: own NMSE [dB] and excess MSE [dB], measured [predicted]\n', res.models{mdl});
    fprintf('    %-16s %8s', 'held out', 'own');
    fprintf(' %18s', res.dictNames{:}); fprintf('\n');
    for j = 1:numel(res.nets)
        own = mean(res.loo.own(j, :, mdl)); den = mean(res.loo.den(j, :, mdl));
        fprintf('    %-16s %8.2f', res.netLabels{j}, 10 * log10(own / den));
        for d = 1:numel(res.dictNames)
            m = 10 * log10(mean(res.loo.meas(j, d, :, mdl)) / own);
            p = 10 * log10(mean(res.loo.pred(j, d, :, mdl)) / own);
            fprintf('   %7.2f [%7.2f]', m, p);
        end
        fprintf('\n');
    end
    fprintf('    max relative difference, EVD-free vs. EVD filters: %.1e\n', ...
            max(reshape(abs(res.loo.meas(:, :, :, mdl) - res.loo.evd(:, :, :, mdl)) ./ ...
                        res.loo.evd(:, :, :, mdl), 1, [])));
    fprintf('  %s few-shot: excess MSE [dB] averaged over the six networks, T0 = %s\n', ...
            res.models{mdl}, mat2str(res.T0));
    for k = 1:numel(res.fewNames)
        fprintf('    %-12s%s\n', res.fewNames{k}, sprintf(' %7.2f', fewExcess(res, k, mdl)));
    end
end
end

function ex = fewExcess(res, k, mdl)
% Excess over the target's own GSP-LMMSE estimator (full training set), in dB,
% averaged over the held-out networks.
J = numel(res.nets); nT0 = numel(res.T0);
ex = zeros(1, nT0);
for t = 1:nT0
    v = zeros(1, J);
    for j = 1:J
        v(j) = mean(reshape(res.few.mse(j, k, t, :, :, mdl), 1, [])) / mean(res.loo.own(j, :, mdl));
    end
    ex(t) = mean(10 * log10(v));
end
end
