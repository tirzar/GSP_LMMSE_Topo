function [est, H] = est_multibasis(fam, S, idx)
%EST_MULTIBASIS  Multi-basis GSP-LMMSE estimator, Theorem 1, eqs. (19)-(20).
%
%   [est, H] = EST_MULTIBASIS(fam, S, idx) uses the graph Fourier bases of the
%   candidates idx (default: all M) and the sample mixture moments S
%   (DAT_MOMENTS), and returns the estimator and the N x numel(idx) responses
%   H = [h_1, ..., h_K] (minimum-norm solution eta = S_D^+ r_D, computed with a
%   negligible ridge 1e-10 tr(S_D)/(KN) for numerical stability).

if nargin < 3 || isempty(idx), idx = 1:fam.M; end
V = cellfun(@(g) g.V, fam.grd(idx), 'UniformOutput', false);
[K, r] = est_spectral_blocks(V, S.CyyMix, S.CxyMix);
SD = cell2mat(K);
n = size(SD, 1);
eta = (SD + 1e-10 * trace(SD) / n * eye(n)) \ r(:);
H = reshape(eta, fam.N, numel(idx));
A = zeros(fam.N);
for i = 1:numel(idx)
    A = A + V{i} * (H(:, i) .* V{i}');
end
key = 'mb';
if numel(idx) < fam.M, key = sprintf('mb_K%d', numel(idx)); end
est = est_affine(key, A, S.muxMix, S.muyMix, numel(idx) * fam.N);   % nominal MN (Table II)
end
