function est = est_single_basis(V, S, key)
%EST_SINGLE_BASIS  Topology-blind graph filter in one fixed basis V.
%
%   est = EST_SINGLE_BASIS(V, S, key) returns the MSE-optimal (for the
%   mixture) estimator x_hat = mux + V diag(h) V' (y - muy), whose response is
%   h = diag(V' Cxy V) ./ diag(V' Cyy V) (Proposition 1 with one basis).
%   Used with V = V_1 (nominal topology, mixture design) and with the basis of
%   the expected Laplacian sum_m p_m L_m (expected-graph filter, cf. GSP over
%   random graphs [Isufi et al. 2017]).

dxy = sum(V .* (S.CxyMix * V), 1)';
dyy = sum(V .* (S.CyyMix * V), 1)';
h = dxy ./ dyy;
est = est_affine(key, V * (h .* V'), S.muxMix, S.muyMix, size(V, 1));
end
