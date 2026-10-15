function [m, lo, hi] = ut_mc_rows(V)
%UT_MC_ROWS  Row-wise Monte Carlo mean and 95% CI [dB], ignoring NaN entries.
%
%   [m, lo, hi] = UT_MC_ROWS(V): V is S x R (points x realizations) in linear
%   scale; m, lo, hi are 1 x S in dB (UT_MC_SUMMARY). Rows without any value
%   give NaN (no marker is drawn).
S = size(V, 1);
m = NaN(1, S); lo = m; hi = m;
for s = 1:S
    v = V(s, :);
    v = v(~isnan(v));
    if isempty(v), continue; end
    [m(s), lo(s), hi(s)] = ut_mc_summary(v, 2);
end
end
