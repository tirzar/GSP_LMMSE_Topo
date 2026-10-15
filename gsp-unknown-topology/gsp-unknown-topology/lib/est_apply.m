function [Xh, mHat] = est_apply(est, Y, mTrue, nSnap)
%EST_APPLY  Apply an estimator to the columns of Y.
%
%   [Xh, mHat] = EST_APPLY(est, Y, mTrue, nSnap)
%     affine  x_hat = A y + b (topology blind)
%     bank    branch of the true topology mTrue (topology-aware benchmark)
%     mm      selective / multiple-model estimator; the columns of Y are
%             grouped into consecutive blocks of nSnap snapshots (default 1)
%             that share the topology, whose log-likelihoods are summed; mHat
%             is the selected topology per column
%     dnn     trained network (EST_DNN)

if nargin < 4 || isempty(nSnap), nSnap = 1; end
mHat = [];
switch est.kind
    case 'affine'
        Xh = est.A * Y + est.b;
    case 'bank'
        Xh = est.A{mTrue} * Y + est.b{mTrue};
    case 'dnn'
        P = est.P;
        Yn = (Y - est.muY) ./ est.sdY;
        out = P.W3 * max(P.W2 * max(P.W1 * Yn + P.b1, 0) + P.b2, 0) + P.b3;
        Xh = out .* est.sdX + est.muX;
    case 'mm'
        [Xh, mHat] = applyMM(est, Y, nSnap);
    otherwise
        error('est_apply:kind', 'Unknown estimator kind "%s".', est.kind);
end
end

function [Xh, mHat] = applyMM(est, Y, nSnap)
M = numel(est.A);
[N, T] = size(Y);
logL = zeros(M, T);
for m = 1:M
    lk = est.lik{m};
    Yc = Y - lk.mu;
    switch lk.type
        case 'diag'
            Z = lk.V' * Yc;
            logL(m, :) = -0.5 * sum(log(lk.d)) - 0.5 * sum(Z .^ 2 ./ lk.d, 1);
        case 'full'
            Z = lk.R' \ Yc;
            logL(m, :) = -sum(log(diag(lk.R))) - 0.5 * sum(Z .^ 2, 1);
    end
end
nB = ceil(T / nSnap);
blk = repelem(1:nB, nSnap);
blk = blk(1:T);
logP = zeros(M, nB);
for m = 1:M
    logP(m, :) = accumarray(blk(:), logL(m, :)', [nB 1])';
end
logP = logP + log(est.p(:));
logP = logP - max(logP, [], 1);
W = exp(logP);
W = W ./ sum(W, 1);
[~, mB] = max(W, [], 1);
mHat = mB(blk);
Xh = zeros(N, T);
if strcmp(est.rule, 'soft')
    Wcol = W(:, blk);
    for m = 1:M
        Xh = Xh + Wcol(m, :) .* (est.A{m} * Y + est.b{m});
    end
else
    for m = 1:M
        idx = (mHat == m);
        if any(idx)
            Xh(:, idx) = est.A{m} * Y(:, idx) + est.b{m};
        end
    end
end
end
