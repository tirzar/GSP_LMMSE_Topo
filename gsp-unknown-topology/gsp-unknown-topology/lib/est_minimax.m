function [est, nu] = est_minimax(fam, S, type, nIter)
%EST_MINIMAX  Minimax-MSE and minimax-regret affine estimators.
%
%   [est, nu] = EST_MINIMAX(fam, S, type) solves, over affine estimators
%   x_hat = A y + b,
%     type = 'mse':     min_{A,b} max_m  MSE_m(A, b)            [Eldar et al.]
%     type = 'regret':  min_{A,b} max_m (MSE_m(A, b) - MSE_m^*) [Eldar & Merhav]
%   where MSE_m is the MSE under topology m (sample conditional moments) and
%   MSE_m^* is that of the topology-aware LMMSE estimator (14). These robust
%   estimators were developed for model uncertainty in linear models; here the
%   uncertainty set is the finite candidate family. Since the objective is
%   convex in (A, b) and linear in nu, the saddle point satisfies
%     A(nu) = Cxy(nu) Cyy(nu)^{-1},   nu* = argmax_{nu in simplex} g(nu),
%   with (A(nu), b(nu)) the LMMSE estimator of the nu-weighted mixture and
%   g(nu) = sum_m nu_m (MSE_m(A(nu), b(nu)) - c_m) concave. nu* is found by
%   mirror ascent (exponentiated gradient) with step halving; the returned
%   estimator is A(nu*). nIter default 300.

if nargin < 4, nIter = 300; end
M = fam.M;
N = fam.N;
c = zeros(M, 1);
Cyx = cellfun(@(C) C', S.Cxy, 'UniformOutput', false);
if strcmp(type, 'regret')                  % MSE of the topology-aware LMMSE
    for m = 1:M                            % (cross-validated ridge, robust for T < N)
        Am = est_ridge(S.Xraw{m} - S.mux(:, m), S.Yraw{m} - S.muy(:, m));
        c(m) = S.trCxx(m) - 2 * sum(sum(Am .* S.Cxy{m})) + sum(sum((Am * S.Cyy{m}) .* Am));
    end
end
nu = ones(M, 1) / M;
[A, b, mse] = designAndEvaluate(nu, S, Cyx, M, N);
g = nu' * (mse - c);
eta = 1 / max(abs(mse - c));
for it = 1:nIter
    nuNew = nu .* exp(eta * (mse - c));
    nuNew = nuNew / sum(nuNew);
    [Anew, bnew, mseNew] = designAndEvaluate(nuNew, S, Cyx, M, N);
    gNew = nuNew' * (mseNew - c);
    if gNew >= g - 1e-12 * abs(g)
        nu = nuNew; A = Anew; b = bnew; mse = mseNew;
        if abs(gNew - g) < 1e-9 * abs(g) && it > 20, g = gNew; break; end
        g = gNew; eta = eta * 1.2;
    else
        eta = eta / 2;
    end
end
est = struct('key', ['minimax_' type], 'kind', 'affine', 'A', A, 'b', b, 'nParams', N ^ 2);
end

function [A, b, mse] = designAndEvaluate(w, S, Cyx, M, N)
    mux = S.mux * w;
    muy = S.muy * w;
    Cxy = zeros(N); Cyy = zeros(N);
    for k = 1:M
        dx = S.mux(:, k) - mux; dy = S.muy(:, k) - muy;
        Cxy = Cxy + w(k) * (S.Cxy{k} + dx * dy');
        Cyy = Cyy + w(k) * (S.Cyy{k} + dy * dy');
    end
    A = Cxy * pinv((Cyy + Cyy') / 2);
    b = mux - A * muy;
    mse = zeros(M, 1);
    for k = 1:M
        e = S.mux(:, k) - A * S.muy(:, k) - b;
        mse(k) = S.trCxx(k) - 2 * sum(sum(A .* Cyx{k}')) + sum(sum((A * S.Cyy{k}) .* A)) + e' * e;
    end
end
