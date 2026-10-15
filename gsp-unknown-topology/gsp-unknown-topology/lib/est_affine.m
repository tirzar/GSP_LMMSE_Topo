function est = est_affine(key, A, mux, muy, nParams)
%EST_AFFINE  Topology-blind affine estimator x_hat = mux + A (y - muy).
%
%   est = EST_AFFINE(key, A, mux, muy, nParams) packs an operator A and the
%   centering means into an estimator struct (kind 'affine').

est = struct('key', key, 'kind', 'affine', 'A', A, 'b', mux - A * muy, ...
             'nParams', nParams);
end
