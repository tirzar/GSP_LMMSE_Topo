function [mse, mseFilt] = xs_test_mse(grd, resp, h, Xte, Yte, dict, alpha, j)
%XS_TEST_MSE  Test MSE of x_hat = mux + V diag(h) V' (y - muy) (eq. (57)).
%
%   mse = XS_TEST_MSE(grd, resp, h, Xte, Yte) uses the EVD (any response h);
%   [mse, mseFilt] = XS_TEST_MSE(..., dict, alpha, j) also returns the MSE of
%   the EVD-free implementation XS_FILTER (identical up to round-off).

Yc = Yte - resp.muy;
Xh = resp.mux + grd.V * (h(:) .* (grd.V' * Yc));
mse = mean(sum((Xh - Xte) .^ 2, 1));
if nargout > 1
    Xf = resp.mux + xs_filter(dict, alpha, grd.L, j, Yc);
    mseFilt = mean(sum((Xf - Xte) .^ 2, 1));
end
end
