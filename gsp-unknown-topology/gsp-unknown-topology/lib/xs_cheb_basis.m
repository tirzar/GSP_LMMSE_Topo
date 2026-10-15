function Psi = xs_cheb_basis(t, Q)
%XS_CHEB_BASIS  Chebyshev polynomials T_0, ..., T_{Q-1} at the points t.
%   Psi(:, q+1) = T_q(t) by the three-term recursion; bounded by 1 on [-1, 1].
t = t(:);
Psi = zeros(numel(t), Q);
Psi(:, 1) = 1;
if Q >= 2, Psi(:, 2) = t; end
for q = 3:Q
    Psi(:, q) = 2 * t .* Psi(:, q - 1) - Psi(:, q - 2);
end
end
