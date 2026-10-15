function ut_seed(s)
%UT_SEED  Seed the random number generators (MATLAB and GNU Octave).
%
%   UT_SEED(s) makes every Monte Carlo realization reproducible: the seed of a
%   realization is a deterministic function of the experiment, the network,
%   and the realization index (UT_SEED_OF).

if exist('OCTAVE_VERSION', 'builtin')
    rand('state', s);  %#ok<RAND>
    randn('state', s); %#ok<RAND>
else
    rng(s, 'twister');
end
end
