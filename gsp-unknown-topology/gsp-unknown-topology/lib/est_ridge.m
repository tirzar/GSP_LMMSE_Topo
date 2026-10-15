function [A, gamma, Cyy] = est_ridge(Xc, Yc, w, nFold)
%EST_RIDGE  LMMSE operator with a cross-validated ridge (diagonal loading).
%
%   [A, gamma, Cyy] = EST_RIDGE(Xc, Yc, w, nFold) returns
%     A = Cxy (Cyy + gamma I)^{-1},  Cxy = sum_t w_t x_t y_t',  Cyy = sum_t w_t y_t y_t'
%   for centered samples (columns of Xc, Yc; weights w, default uniform) with
%   gamma = g * tr(Cyy)/N, g chosen from logspace(-8, 0, 17) by nFold-fold
%   cross-validation of the estimation MSE (default 5 folds). Cyy is returned
%   without loading.

[N, T] = size(Yc);
if nargin < 3 || isempty(w), w = ones(1, T) / T; end
if nargin < 4 || isempty(nFold), nFold = 5; end
w = w(:)' / sum(w);
fold = mod(randperm(T), nFold) + 1;
grid = logspace(-8, 0, 17);
err = zeros(1, numel(grid));
for f = 1:nFold
    tr = fold ~= f; te = ~tr;
    wtr = w(tr) / sum(w(tr));
    Cxy = (Xc(:, tr) .* wtr) * Yc(:, tr)';
    Cyy = (Yc(:, tr) .* wtr) * Yc(:, tr)';
    [U, Dg] = eig((Cyy + Cyy') / 2);
    d = max(real(diag(Dg)), 0);
    g0 = trace(Cyy) / N;
    CU = Cxy * U;
    Z = U' * Yc(:, te);
    for k = 1:numel(grid)
        E = Xc(:, te) - CU * (Z ./ (d + grid(k) * g0));
        err(k) = err(k) + sum(w(te) .* sum(E .^ 2, 1));
    end
end
[~, k] = min(err);
Cxy = (Xc .* w) * Yc';
Cyy = (Yc .* w) * Yc';
Cyy = (Cyy + Cyy') / 2;
gamma = grid(k) * trace(Cyy) / N;
A = Cxy / (Cyy + gamma * eye(N));
end
