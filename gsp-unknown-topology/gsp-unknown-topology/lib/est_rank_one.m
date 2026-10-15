function est = est_rank_one(fam, S, nIter)
%EST_RANK_ONE  Rank-one multi-basis estimator, Remark 1 (r = 1), learned weights.
%
%   est = EST_RANK_ONE(fam, S) fits H = h b' (h: common response, b: topology
%   weights normalized to 1'b = 1) by alternating least squares: for fixed b,
%   h is given by Proposition 1; for fixed h, b solves an M x M least-squares
%   problem. Initialized at b = p; nIter iterations (default 50).

if nargin < 3, nIter = 50; end
V = cellfun(@(g) g.V, fam.grd, 'UniformOutput', false);
[K, r] = est_spectral_blocks(V, S.CyyMix, S.CxyMix);
M = fam.M;
N = fam.N;
b = fam.p(:);
for it = 1:nIter
    h = sharedResponse(K, r, b, N);
    G = zeros(M); g = zeros(M, 1);
    for m = 1:M
        Kh = cellfun(@(Kml) Kml * h, K(m, :), 'UniformOutput', false);
        G(m, :) = h' * cell2mat(Kh);
        g(m) = h' * r(:, m);
    end
    bNew = (G + 1e-10 * trace(G) / M * eye(M)) \ g;
    bNew = bNew / sum(bNew);
    done = norm(bNew - b) < 1e-8 * norm(b);
    b = bNew;
    if done, break; end
end
h = sharedResponse(K, r, b, N);
A = zeros(N);
for m = 1:M
    A = A + b(m) * V{m} * (h .* V{m}');
end
est = est_affine('r1', A, S.muxMix, S.muyMix, N + M - 1);
end

function h = sharedResponse(K, r, w, N)
Ssr = zeros(N);
for m = 1:numel(w)
    for l = 1:numel(w)
        Ssr = Ssr + w(m) * w(l) * K{m, l};
    end
end
h = (Ssr + 1e-10 * trace(Ssr) / N * eye(N)) \ (r * w(:));
end
