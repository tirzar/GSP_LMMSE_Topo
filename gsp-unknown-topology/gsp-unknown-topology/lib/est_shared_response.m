function est = est_shared_response(fam, S, w, key)
%EST_SHARED_RESPONSE  Shared-response estimator, Proposition 1, eqs. (29)-(30).
%
%   est = EST_SHARED_RESPONSE(fam, S, w) designs the common response h for the
%   fixed weights w (M x 1), x_hat = mux + sum_m w_m V_m diag(h) V_m' (y - muy).
%   w = fam.p gives prior-weighted fusion; w = e_m0 gives the topology-blind
%   single-basis estimator of candidate m0.

if nargin < 4, key = 'sr'; end
idx = find(w ~= 0);
V = cellfun(@(g) g.V, fam.grd(idx), 'UniformOutput', false);
[K, r] = est_spectral_blocks(V, S.CyyMix, S.CxyMix);
w = w(idx);
Ssr = zeros(fam.N);
for m = 1:numel(idx)
    for l = 1:numel(idx)
        Ssr = Ssr + w(m) * w(l) * K{m, l};
    end
end
h = (Ssr + 1e-10 * trace(Ssr) / fam.N * eye(fam.N)) \ (r * w(:));
A = zeros(fam.N);
for m = 1:numel(idx)
    A = A + w(m) * V{m} * (h .* V{m}');
end
est = est_affine(key, A, S.muxMix, S.muyMix, fam.N);
end
