function net = xs_network_data(cfg, caseName, model, T, Ttest, seed)
%XS_NETWORK_DATA  One network with a known topology and its data (Section VI-F).
%
%   net = XS_NETWORK_DATA(cfg, caseName, model, T, Ttest, seed) builds the graph
%   (Laplacian QA enforced) and draws T training and Ttest test samples of
%   data model 'M1' or 'M2': fields grd, X, G, W (unit noise), Xte, Gte, Wte.
%   Measurements at noise level sigma: y = G + sigma W.

mpc = net_load_case(caseName);
grd = net_build_grid(mpc);
[ok, qa] = net_check_laplacian(grd, mpc);
if ~ok
    disp(qa);
    error('xs_network_data:qa', 'Laplacian QA failed for %s.', caseName);
end
ut_seed(seed);
switch upper(model)
    case 'M1'
        X = dat_prior_states(grd, cfg.beta, T + Ttest);
        G = dat_ac_injection(grd.Ybus, X);
    case 'M2'
        [X, G] = dat_pf_states(mpc, grd, T + Ttest, cfg.pf);
end
W = randn(grd.N, T + Ttest);
net = struct('name', caseName, 'grd', grd, 'X', X(:, 1:T), 'G', G(:, 1:T), 'W', W(:, 1:T), ...
             'Xte', X(:, T+1:end), 'Gte', G(:, T+1:end), 'Wte', W(:, T+1:end));
end
