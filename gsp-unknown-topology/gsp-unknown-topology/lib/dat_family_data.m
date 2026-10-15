function D = dat_family_data(fam, model, Ttrain, Ttest, cfg)
%DAT_FAMILY_DATA  Training and test data for every candidate topology.
%
%   D = DAT_FAMILY_DATA(fam, model, Ttrain, Ttest, cfg) draws, for each
%   candidate m of the family fam (NET_TOPOLOGY_FAMILY):
%     D.Xtr{m}, D.Gtr{m}, D.Wtr{m}   N x Ttrain states, noise-free injections
%                                    g(L_m, x) (eq. (70)), unit white noise
%     D.Xte{m}, D.Gte{m}, D.Wte{m}   the same for Ttest test samples
%   model = 'M1' (smooth prior, cfg.beta) or 'M2' (AC power-flow states,
%   cfg.pf). Measurements at noise level sigma are y = G + sigma * W
%   (common random numbers across noise levels).

D.M = fam.M;
D.N = fam.N;
D.model = model;
D.Xtr = cell(1, fam.M); D.Gtr = D.Xtr; D.Wtr = D.Xtr;
D.Xte = D.Xtr; D.Gte = D.Xtr; D.Wte = D.Xtr;
for m = 1:fam.M
    g = fam.grd{m};
    switch upper(model)
        case 'M1'
            X = dat_prior_states(g, cfg.beta, Ttrain + Ttest);
            G = dat_ac_injection(g.Ybus, X);
        case 'M2'
            [X, G] = dat_pf_states(fam.mpc{m}, g, Ttrain + Ttest, cfg.pf);
        otherwise
            error('dat_family_data:model', 'Unknown data model "%s".', model);
    end
    W = randn(g.N, Ttrain + Ttest);
    D.Xtr{m} = X(:, 1:Ttrain);      D.Xte{m} = X(:, Ttrain+1:end);
    D.Gtr{m} = G(:, 1:Ttrain);      D.Gte{m} = G(:, Ttrain+1:end);
    D.Wtr{m} = W(:, 1:Ttrain);      D.Wte{m} = W(:, Ttrain+1:end);
end
end
