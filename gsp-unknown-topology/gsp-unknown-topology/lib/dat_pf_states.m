function [X, G] = dat_pf_states(mpc, grd, T, opts)
%DAT_PF_STATES  Data model M2: phase angles from AC power-flow solutions.
%
%   [X, G] = DAT_PF_STATES(mpc, grd, T, opts) draws T operating points:
%     loads       P_d,n = kappa * P_d,n^base * (1 + opts.sigLoad * eps_n),
%                 kappa ~ U[opts.kappa], eps_n ~ N(0, 1)
%     generation  shared by the in-service generators in proportion to their
%                 base dispatch; the slack bus balances the losses
%     solution    Newton-Raphson on the active power of the non-slack buses
%                 (|v| = 1 p.u.), initialized by the DC power flow
%   X = mean-centered angles (angles are defined up to a common shift), and
%   G = g(L, x), the resulting AC injections at all buses (eq. (70)).
%   Unlike M1, the state distribution is not defined through L, and the
%   injections depend on the topology only through the losses.

if ~isfield(opts, 'sigLoad'), opts.sigLoad = 0.1; end
if ~isfield(opts, 'kappa'),   opts.kappa = [0.8 1.2]; end
N = grd.N;
Y = grd.Ybus;
e2i = zeros(max(mpc.bus(:, 1)), 1);
e2i(mpc.bus(:, 1)) = (1:N)';
ns = setdiff(1:N, grd.ref);
Pd0 = mpc.bus(:, 3) / mpc.baseMVA;
on = mpc.gen(:, 8) > 0;
pg0 = max(mpc.gen(on, 2), 0);
if sum(pg0) <= 0, pg0 = mpc.gen(on, 9); end
share = accumarray(e2i(mpc.gen(on, 1)), pg0, [N 1]) / sum(pg0);
Lns = grd.L(ns, ns);

X = zeros(N, T);
G = zeros(N, T);
t = 1;
rejected = 0;
while t <= T
    kappa = opts.kappa(1) + diff(opts.kappa) * rand;
    Pd = kappa * Pd0 .* (1 + opts.sigLoad * randn(N, 1));
    Pspec = share * sum(Pd) - Pd;
    th = zeros(N, 1);
    th(ns) = Lns \ Pspec(ns);
    ok = false;
    for it = 1:25
        V = exp(1j * th);
        I = Y * V;
        F = real(V .* conj(I)) - Pspec;
        if max(abs(F(ns))) < 1e-9, ok = true; break; end
        J = real(1j * diag(V) * conj(diag(I) - Y * diag(V)));
        th(ns) = th(ns) - J(ns, ns) \ F(ns);
    end
    if ~ok || any(~isfinite(th))
        rejected = rejected + 1;
        if rejected > 10 * T
            error('dat_pf_states:diverge', 'AC power flow does not converge.');
        end
        continue;
    end
    X(:, t) = th - mean(th);
    G(:, t) = real(V .* conj(I));
    t = t + 1;
end
end
