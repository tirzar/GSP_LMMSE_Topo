function P = dat_ac_injection(Ybus, Theta, Vm)
%DAT_AC_INJECTION  AC active-power injections, eq. (70).
%
%   P = DAT_AC_INJECTION(Ybus, Theta) returns, column by column,
%     [g(L, x)]_n = sum_k |v_n||v_k| (G_nk cos(x_n - x_k) + B_nk sin(x_n - x_k)),
%   with Ybus = G + jB, Theta = N x T phase angles [rad], and |v| = 1 p.u.

if nargin < 3 || isempty(Vm)
    Vm = ones(size(Theta, 1), 1);
end
Vc = Vm(:) .* exp(1j * Theta);
P = real(Vc .* conj(Ybus * Vc));
end
