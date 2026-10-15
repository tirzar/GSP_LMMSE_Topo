function [mDB, loDB, hiDB] = ut_mc_summary(V, dim)
%UT_MC_SUMMARY  Monte Carlo mean [dB] with a simple 95% confidence interval.
%
%   [mDB, loDB, hiDB] = UT_MC_SUMMARY(V, dim) averages the linear-scale values
%   V (e.g., NMSE) over the Monte Carlo dimension dim (default: last) and
%   returns mDB = 10 log10(mean(V)). The interval is mDB -/+ 1.96 standard
%   errors of the per-realization values in dB, i.e.,
%       loDB, hiDB = mDB -/+ 1.96 * std(10 log10 V) / sqrt(R),
%   which is symmetric in dB and never degenerate.

if nargin < 2, dim = ndims(V); end
R = size(V, dim);
mDB = 10 * log10(mean(V, dim));
if R > 1
    half = 1.96 * std(10 * log10(V), 0, dim) / sqrt(R);
else
    half = zeros(size(mDB));
end
loDB = mDB - half;
hiDB = mDB + half;
end
