function dict = xs_dictionary(kind, varargin)
%XS_DICTIONARY  Fixed spectral dictionary {phi_q(lambda)} of Section IV, eq. (45).
%
%   dict = XS_DICTIONARY('chebyshev', Q, s, mu)
%     phi_q(lambda) = T_q(lambda / s_j - 1), q = 0, ..., Q-1 (Chebyshev polynomials
%     of the first kind, i.e., a polynomial graph filter in L/s_j - I, cf. (4)),
%     with the spectral scale s_j of network j (s scalar: common absolute axis;
%     s vector: one scale per network, e.g., s_j = lambda_max,j / 2 for the
%     normalized axis of ChebNet [Defferrard et al.] or s_j proportional to N_j
%     for the graphon-normalized axis lambda/N [Ruiz et al.]), and the stability
%     penalty mu * sum_q q^4 alpha_q^2.
%   dict = XS_DICTIONARY('rational', poles, mu)
%     phi_q(lambda) = 1/(lambda + c_q) with fixed poles c_q > 0 and the ridge
%     penalty mu * sum_q alpha_q^2.
%
%   Fields: kind, Q, s or poles, pen (1 x Q), Phi(lambda, j) (numel(lambda) x Q).

switch lower(kind)
    case 'chebyshev'
        Q = varargin{1}; s = varargin{2}; mu = varargin{3};
        dict = struct('kind', 'chebyshev', 'Q', Q, 's', s(:)', 'poles', [], ...
            'pen', mu * (0:Q-1) .^ 4);
        dict.Phi = @(lam, j) xs_cheb_basis(lam(:) / s(min(j, numel(s))) - 1, Q);
    case 'rational'
        c = varargin{1}(:)'; mu = varargin{2};
        dict = struct('kind', 'rational', 'Q', numel(c), 's', [], 'poles', c, ...
            'pen', mu * ones(1, numel(c)));
        dict.Phi = @(lam, j) 1 ./ (lam(:) + c);
    otherwise
        error('xs_dictionary:kind', 'Unknown dictionary "%s".', kind);
end
end
