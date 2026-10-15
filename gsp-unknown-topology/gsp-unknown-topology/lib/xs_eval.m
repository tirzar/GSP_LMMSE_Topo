function mse = xs_eval(te, V, mux, muy, h)
%XS_EVAL  Exact test MSE of a graph filter from the test-set GFT moments.
%
%   mse = XS_EVAL(te, V, mux, muy, h) returns the average of
%   ||mux + V diag(h) V' (y - muy) - x||^2 over the test set whose moments te
%   (XS_MOMENTS of the test data) were computed in the basis V. Since the error
%   decouples across graph frequencies,
%     MSE = sum_n (dx_n - h_n dy_n)^2 + h_n^2 Dyy_n - 2 h_n Dxy_n + Dxx_n,
%   with dx = V'(mux - mux_te), dy = V'(muy - muy_te) and the test moments
%   Dxx, Dxy, Dyy: exact, and O(N) per response once the moments are known.

h = h(:);
dx = V' * (mux - te.mux);
dy = V' * (muy - te.muy);
mse = sum((dx - h .* dy) .^ 2 + h .^ 2 .* te.dyy - 2 * h .* te.dxy + te.dxx);
end
