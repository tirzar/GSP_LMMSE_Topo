function [famSeed, dataSeed] = ut_seed_of(cfg, netIdx, r, expId, variant)
%UT_SEED_OF  Seeds of Monte Carlo realization r of network netIdx.
%
%   [famSeed, dataSeed] = UT_SEED_OF(cfg, netIdx, r, expId, variant)
%   famSeed depends only on (netIdx, r, variant), so that Sections VI-B to VI-E
%   use the same candidate families; dataSeed also depends on the experiment
%   expId (2..5 for Sections VI-B to VI-E). variant (default 0) distinguishes the perturbation
%   sizes of Section VI-C.

if nargin < 5 || isempty(variant), variant = 0; end
famSeed  = cfg.seed + 1000 * netIdx + r + 100000 * variant;
dataSeed = famSeed + 5000000 + 1000000 * expId;
end
