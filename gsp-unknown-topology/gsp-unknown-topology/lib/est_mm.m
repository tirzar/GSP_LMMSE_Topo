function est = est_mm(fam, S, branch, rule)
%EST_MM  Topology-selective and multiple-model (MM) estimators (Sec. III-D).
%
%   est = EST_MM(fam, S, branch, rule) combines topology-conditioned affine
%   branches x_hat_m(y) = mux_m + A_m (y - muy_m) using approximate posterior
%   topology probabilities pi_m(y) ~ p_m p_hat(y | m):
%     branch = 'gsp'   GSP-LMMSE branches (33) and the graph-stationary Gaussian
%                      surrogate likelihood N(muy_m, V_m diag(dyy_m) V_m')
%                      (diagonal in the GFT domain, O(N^2) per candidate)
%     branch = 'full'  sample LMMSE branches (14) and the full Gaussian
%                      likelihood N(muy_m, Cyy_m + gamma_m I), both with the
%                      cross-validated loading gamma_m of EST_RIDGE: the
%                      classical multiple-model (Magill) estimator
%     rule = 'hard'    x_hat = x_hat_{m_hat}(y), m_hat = argmax pi_m (eq. (38))
%     rule = 'soft'    x_hat = sum_m pi_m(y) x_hat_m(y), eq. (40), which
%                      approximates the topology-blind MMSE estimator (11)
%   With P snapshots that share the topology, the log-likelihoods are summed
%   before the decision (EST_APPLY).

M = fam.M;
N = fam.N;
A = cell(1, M); b = cell(1, M); lik = cell(1, M);
for m = 1:M
    switch branch
        case 'gsp'
            V = fam.grd{m}.V;
            A{m} = V * ((S.dxy(:, m) ./ S.dyy(:, m)) .* V');
            lik{m} = struct('type', 'diag', 'V', V, 'd', S.dyy(:, m), 'mu', S.muy(:, m));
        case 'full'
            [A{m}, gam, Cyy] = est_ridge(S.Xraw{m} - S.mux(:, m), S.Yraw{m} - S.muy(:, m));
            R = chol(Cyy + max(gam, 1e-10 * trace(Cyy) / N) * eye(N));
            lik{m} = struct('type', 'full', 'R', R, 'mu', S.muy(:, m));
    end
    b{m} = S.mux(:, m) - A{m} * S.muy(:, m);
end
if strcmp(branch, 'gsp'), np = M * N; else, np = M * N ^ 2; end
est = struct('key', ['mm_' branch '_' rule], 'kind', 'mm', 'A', {A}, 'b', {b}, ...
             'lik', {lik}, 'p', S.p, 'rule', rule, 'nParams', np);
end
