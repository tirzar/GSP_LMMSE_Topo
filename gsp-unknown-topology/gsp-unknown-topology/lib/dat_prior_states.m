function X = dat_prior_states(grd, beta, T)
%DAT_PRIOR_STATES  Data model M1: smooth Gaussian graph-signal prior.
%
%   X = DAT_PRIOR_STATES(grd, beta, T) draws T samples x ~ N(0, beta * L^+),
%   i.e., x~_1 = 0 and x~_n ~ N(0, beta/lambda_n) for n >= 2 [Kroizer et al.].

N = grd.N;
X = grd.V(:, 2:N) * (sqrt(beta ./ grd.lambda(2:N)) .* randn(N - 1, T));
end
