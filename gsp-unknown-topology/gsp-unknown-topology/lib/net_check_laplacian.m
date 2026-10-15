function [ok, qa] = net_check_laplacian(grd, mpc)
%NET_CHECK_LAPLACIAN  Validity checks of a graph Laplacian and its EVD.
%
%   [ok, qa] = NET_CHECK_LAPLACIAN(grd, mpc) verifies that grd.L is symmetric,
%   has zero row sums and nonpositive off-diagonal entries, is PSD with a simple
%   zero eigenvalue (connected graph), has distinct eigenvalues (as assumed in
%   Sec. II-B of the paper), that the stored EVD is accurate, and that the edges
%   of L are the in-service branches of mpc. qa also reports the relative
%   distance between L and the finite-difference Jacobian dP/dtheta at theta = 0.

L = grd.L;
N = size(L, 1);
tol = 1e-9 * max(abs(L(:)));
off = L - diag(diag(L));
lam = sort(eig((L + L') / 2));

qa.symmetric = max(max(abs(L - L'))) <= tol;
qa.rowSum = max(abs(sum(L, 2))) <= tol;
qa.offDiag = all(off(:) <= tol);
qa.psd = lam(1) >= -tol;
qa.connected = lam(2) > tol;
qa.distinct = min(diff(lam(2:end))) > 1e-9 * lam(end);
qa.evd = norm(grd.V * diag(grd.lambda) * grd.V' - L, 'fro') <= 1e-8 * norm(L, 'fro');

e2i = zeros(max(mpc.bus(:, 1)), 1);
e2i(mpc.bus(:, 1)) = (1:N)';
br = mpc.branch(mpc.branch(:, 11) > 0, :);
A = full(sparse(e2i(br(:, 1)), e2i(br(:, 2)), 1, N, N));
A = (A + A') > 0;
A(1:N+1:end) = false;
qa.support = isequal(A, abs(off) > tol);

E = 1e-6 * eye(N);
J = (dat_ac_injection(grd.Ybus, E) - dat_ac_injection(grd.Ybus, -E)) / 2e-6;
qa.jacobianErr = norm(L - J, 'fro') / norm(J, 'fro');

ok = qa.symmetric && qa.rowSum && qa.offDiag && qa.psd && qa.connected && ...
     qa.distinct && qa.evd && qa.support;
end
