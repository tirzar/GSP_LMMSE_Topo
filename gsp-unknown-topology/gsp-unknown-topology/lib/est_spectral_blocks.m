function [K, r] = est_spectral_blocks(Vlist, Cyy, Cxy)
%EST_SPECTRAL_BLOCKS  Blocks of S_D and r_D in eqs. (17)-(18).
%
%   [K, r] = EST_SPECTRAL_BLOCKS(Vlist, Cyy, Cxy) returns, for the bases
%   Vlist{1..M}, the M x M cell K with K{m,l} = (V_m' Cyy V_l) .* (V_m' V_l)
%   and the N x M matrix r with r(:,m) = diag(V_m' Cxy V_m).

M = numel(Vlist);
N = size(Vlist{1}, 1);
K = cell(M, M);
r = zeros(N, M);
U = cellfun(@(V) Cyy * V, Vlist, 'UniformOutput', false);
for m = 1:M
    for l = m:M
        K{m, l} = (Vlist{m}' * U{l}) .* (Vlist{m}' * Vlist{l});
        K{l, m} = K{m, l}';
    end
    r(:, m) = sum(Vlist{m} .* (Cxy * Vlist{m}), 1)';
end
end
