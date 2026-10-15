function est = est_dnn(S, opt)
%EST_DNN  Data-driven Bayesian state estimator (fully connected DNN).
%
%   est = EST_DNN(S, opt) trains a multilayer perceptron y -> x on the pooled
%   labeled training data of all candidate topologies, without topology
%   labels, following the deep-learning Bayesian state estimator of Mestav,
%   Luengo-Rozas, and Tong (IEEE Trans. Power Syst., 2019). The network
%   approximates E[x | y] by minimizing the empirical MSE.
%   Architecture: N -> opt.hidden(1) -> opt.hidden(2) -> N, ReLU, standardized
%   inputs and outputs. Training: Adam (opt.lr), mini-batches (opt.batch),
%   weight decay (opt.wd), and early stopping on a validation split
%   (opt.valFrac, opt.patience, opt.maxEpochs). No toolbox is required.

if ~isfield(opt, 'hidden'),    opt.hidden = [256 256]; end
if ~isfield(opt, 'lr'),        opt.lr = 1e-3; end
if ~isfield(opt, 'batch'),     opt.batch = 128; end
if ~isfield(opt, 'wd'),        opt.wd = 1e-5; end
if ~isfield(opt, 'valFrac'),   opt.valFrac = 0.15; end
if ~isfield(opt, 'patience'),  opt.patience = 15; end
if ~isfield(opt, 'maxEpochs'), opt.maxEpochs = 200; end

X = cell2mat(S.Xraw);
Y = cell2mat(S.Yraw);
[N, n] = size(Y);
muY = mean(Y, 2); sdY = std(Y, 0, 2) + 1e-12;
muX = mean(X, 2); sdX = std(X, 0, 2) + 1e-12;
Yn = (Y - muY) ./ sdY;
Xn = (X - muX) ./ sdX;
perm = randperm(n);
nVal = max(1, round(opt.valFrac * n));
iVal = perm(1:nVal); iTr = perm(nVal+1:end);
H1 = opt.hidden(1); H2 = opt.hidden(2);
P.W1 = randn(H1, N) * sqrt(2 / N);   P.b1 = zeros(H1, 1);
P.W2 = randn(H2, H1) * sqrt(2 / H1); P.b2 = zeros(H2, 1);
P.W3 = randn(N, H2) * sqrt(1 / H2);  P.b3 = zeros(N, 1);
names = fieldnames(P);
for k = 1:numel(names)
    Mo.(names{k}) = zeros(size(P.(names{k})));
    Ve.(names{k}) = zeros(size(P.(names{k})));
end
beta1 = 0.9; beta2 = 0.999; step = 0;
best = inf; bestP = P; wait = 0;
nTr = numel(iTr);
for epoch = 1:opt.maxEpochs
    ord = iTr(randperm(nTr));
    for s = 1:opt.batch:nTr
        idx = ord(s:min(s + opt.batch - 1, nTr));
        B = numel(idx);
        y = Yn(:, idx); x = Xn(:, idx);
        Z1 = P.W1 * y + P.b1;   A1 = max(Z1, 0);
        Z2 = P.W2 * A1 + P.b2;  A2 = max(Z2, 0);
        out = P.W3 * A2 + P.b3;
        dO = 2 * (out - x) / (B * N);
        G.W3 = dO * A2' + opt.wd * P.W3;  G.b3 = sum(dO, 2);
        dA2 = (P.W3' * dO) .* (Z2 > 0);
        G.W2 = dA2 * A1' + opt.wd * P.W2; G.b2 = sum(dA2, 2);
        dA1 = (P.W2' * dA2) .* (Z1 > 0);
        G.W1 = dA1 * y' + opt.wd * P.W1;  G.b1 = sum(dA1, 2);
        step = step + 1;
        for k = 1:numel(names)
            f = names{k};
            Mo.(f) = beta1 * Mo.(f) + (1 - beta1) * G.(f);
            Ve.(f) = beta2 * Ve.(f) + (1 - beta2) * G.(f) .^ 2;
            mh = Mo.(f) / (1 - beta1 ^ step);
            vh = Ve.(f) / (1 - beta2 ^ step);
            P.(f) = P.(f) - opt.lr * mh ./ (sqrt(vh) + 1e-8);
        end
    end
    vl = mean(sum((forward(P, Yn(:, iVal)) - Xn(:, iVal)) .^ 2, 1));
    if vl < best * (1 - 1e-4)
        best = vl; bestP = P; wait = 0;
    else
        wait = wait + 1;
        if wait >= opt.patience, break; end
    end
end
est = struct('key', 'dnn', 'kind', 'dnn', 'P', bestP, 'muY', muY, 'sdY', sdY, ...
             'muX', muX, 'sdX', sdX, 'epochs', epoch, ...
             'nParams', H1 * (N + 1) + H2 * (H1 + 1) + N * (H2 + 1));
end

function out = forward(P, y)
out = P.W3 * max(P.W2 * max(P.W1 * y + P.b1, 0) + P.b2, 0) + P.b3;
end
