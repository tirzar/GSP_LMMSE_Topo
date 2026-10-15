function Z = xs_filter(dict, alpha, L, j, Y)
%XS_FILTER  EVD-free application of h(L; alpha) to the columns of Y (Section IV-B).
%
%   Chebyshev: three-term recursion u_q = 2 (L/s_j - I) u_{q-1} - u_{q-2}, i.e.,
%   only sparse matrix-vector products (the output at a node depends on its
%   (Q-1)-hop neighborhood) [Shuman et al. 2018].
%   Rational: sum_q alpha_q (L + c_q I)^{-1} Y by sparse linear solves (or ARMA
%   graph filters [Isufi et al. 2017]).
%   The DC mode is projected out in both cases (Z <- Z - mean(Z)).

N = size(L, 1);
Ls = sparse(L);
switch dict.kind
    case 'chebyshev'
        Lb = Ls / dict.s(min(j, numel(dict.s))) - speye(N);
        u0 = Y; Z = alpha(1) * u0;
        if dict.Q >= 2, u1 = Lb * Y; Z = Z + alpha(2) * u1; end
        for q = 3:dict.Q
            u2 = 2 * (Lb * u1) - u0;
            Z = Z + alpha(q) * u2;
            u0 = u1; u1 = u2;
        end
    case 'rational'
        Z = zeros(size(Y));
        for q = 1:dict.Q
            Z = Z + alpha(q) * ((Ls + dict.poles(q) * speye(N)) \ Y);
        end
end
Z = Z - mean(Z, 1);
end
