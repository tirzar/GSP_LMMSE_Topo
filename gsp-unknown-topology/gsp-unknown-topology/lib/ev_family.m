function R = ev_family(ests, D, sigma2, p, nSnap)
%EV_FAMILY  Test-set performance of a list of estimators.
%
%   R = EV_FAMILY(ests, D, sigma2, p, nSnap) applies each estimator in the cell
%   ests to the test data of every candidate topology, y = g(L_m, x) + sigma w,
%   with blocks of nSnap snapshots sharing the topology for the selective and
%   multiple-model estimators (default 1).
%     R.mse(e, m)   test MSE under topology m
%     R.nmse(e)     NMSE = sum_m p_m MSE_m / sum_m p_m E_m||x - mu_x||^2, with
%                   mu_x the mixture mean (linear scale)
%     R.worst(e)    worst-topology NMSE max_m MSE_m / E_m||x - mu_x||^2
%     R.id(e)       topology-identification rate sum_m p_m Pr(m_hat = m | m)
%                   (selective and multiple-model estimators; NaN otherwise)

if nargin < 5 || isempty(nSnap), nSnap = 1; end
M = numel(D.Xte);
E = numel(ests);
p = p(:);
muMix = zeros(size(D.Xte{1}, 1), 1);
for m = 1:M, muMix = muMix + p(m) * mean(D.Xte{m}, 2); end
den = zeros(1, M);
for m = 1:M
    den(m) = mean(sum((D.Xte{m} - muMix) .^ 2, 1));
end
mse = zeros(E, M);
id = NaN(E, M);
for m = 1:M
    Y = D.Gte{m} + sqrt(sigma2) * D.Wte{m};
    for e = 1:E
        [Xh, mHat] = est_apply(ests{e}, Y, m, nSnap);
        mse(e, m) = mean(sum((Xh - D.Xte{m}) .^ 2, 1));
        if ~isempty(mHat), id(e, m) = mean(mHat == m); end
    end
end
R.keys = cellfun(@(s) s.key, ests, 'UniformOutput', false);
R.mse = mse;
R.den = den;
R.nmse = (mse * p) / (den * p);
R.worst = max(mse ./ den, [], 2);
R.id = id * p;
R.nParams = cellfun(@(s) s.nParams, ests);
end
