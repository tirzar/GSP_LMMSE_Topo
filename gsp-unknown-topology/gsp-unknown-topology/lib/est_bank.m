function est = est_bank(fam, S, type)
%EST_BANK  Topology-aware benchmarks x_hat = mux_m + A_m (y - muy_m).
%
%   est = EST_BANK(fam, S, type) returns one affine estimator per candidate
%   topology; at test time, the branch of the true topology is applied.
%     'gsp'           GSP-LMMSE estimator (33)-(35): A_m = V_m diag(dxy_m ./ dyy_m) V_m'
%     'lmmse_ridge'   LMMSE estimator (14) with cross-validated diagonal
%                     loading (EST_RIDGE)

M = fam.M;
N = fam.N;
A = cell(1, M);
switch type
    case 'gsp'
        for m = 1:M
            V = fam.grd{m}.V;
            A{m} = V * ((S.dxy(:, m) ./ S.dyy(:, m)) .* V');
        end
        np = N;
    case 'lmmse_ridge'
        for m = 1:M, A{m} = est_ridge(S.Xraw{m} - S.mux(:, m), S.Yraw{m} - S.muy(:, m)); end
        np = N ^ 2;
    otherwise
        error('est_bank:type', 'Unknown type "%s".', type);
end
b = cell(1, M);
for m = 1:M, b{m} = S.mux(:, m) - A{m} * S.muy(:, m); end
est = struct('key', ['aware_' type], 'kind', 'bank', 'A', {A}, 'b', {b}, 'nParams', np);
end
