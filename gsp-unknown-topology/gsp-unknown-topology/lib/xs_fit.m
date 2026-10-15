function alpha = xs_fit(resp, lam, idx, dict)
%XS_FIT  Cross-network spectral coefficients, Theorem 3, (56) and (68).
%
%   alpha = XS_FIT(resp, lam, idx, dict) minimizes over alpha the weighted
%   curve fit over the networks idx (each a single known graph, Remark 4)
%     sum_j (1/N_j) sum_{n >= 2} dyy_j(n) (h(lambda_jn; alpha) - f_j(n))^2
%       + sum_q pen_q alpha_q^2,
%   which, by completing the square, is the average per-node empirical MSE of
%   the cross-network estimator (67) up to a constant, plus the penalty.
%   resp{j} (XS_MOMENTS) holds f_j and dyy_j; lam{j} the graph frequencies; the
%   DC mode (lambda = 0) is excluded. dict: XS_DICTIONARY.

A = zeros(0, dict.Q);
b = zeros(0, 1);
for j = idx(:)'
    N = numel(lam{j});
    k = 2:N;
    w = sqrt(resp{j}.dyy(k) / N);
    A = [A; w .* dict.Phi(lam{j}(k), j)];  %#ok<AGROW>
    b = [b; w .* resp{j}.f(k)];             %#ok<AGROW>
end
A = [A; diag(sqrt(dict.pen))];
b = [b; zeros(dict.Q, 1)];
alpha = A \ b;
end
