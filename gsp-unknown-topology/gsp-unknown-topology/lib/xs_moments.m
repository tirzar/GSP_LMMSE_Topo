function resp = xs_moments(grd, X, Y)
%XS_MOMENTS  Graph-frequency moments of one network (known topology), eq. (55).
%
%   resp = XS_MOMENTS(grd, X, Y) from paired samples (columns of X and Y):
%     mux, muy        sample means (the affine offset)
%     dxx, dxy, dyy   diag(V' Cxx V), diag(V' Cxy V), diag(V' Cyy V)
%     f = dxy ./ dyy  GSP-LMMSE response (f(1) = 0 for the DC mode)
%     mseGSP          its in-sample MSE, trCxx
%   In-sample, MSE(h) = mseGSP + sum_n dyy(n) (h(n) - f(n))^2 for any response
%   h (completing the square; the identity behind the predicted excess MSE).

T = size(X, 2);
mux = mean(X, 2); muy = mean(Y, 2);
Xt = grd.V' * (X - mux);
Yt = grd.V' * (Y - muy);
dxx = sum(Xt .^ 2, 2) / T;
dxy = sum(Xt .* Yt, 2) / T;
dyy = sum(Yt .^ 2, 2) / T;
f = dxy ./ dyy;
f(1) = 0;
resp = struct('mux', mux, 'muy', muy, 'dxx', dxx, 'dxy', dxy, 'dyy', dyy, 'f', f, ...
              'mseGSP', sum(dxx) - sum(dxy(2:end) .^ 2 ./ dyy(2:end)), 'trCxx', sum(dxx), 'T', T);
end
