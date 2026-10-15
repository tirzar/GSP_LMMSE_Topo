function ut_print_summary(res, snrDB, sIdx)
%UT_PRINT_SUMMARY  Console table: mean NMSE and worst-topology NMSE [dB].
%
%   UT_PRINT_SUMMARY(res, snrDB) prints, for each network of res and the SNR
%   snrDB, the Monte Carlo mean (and 95% CI half-width) of the NMSE, the
%   worst-topology NMSE, and the identification rate of every method.
%   UT_PRINT_SUMMARY(res, [], sIdx) uses the index sIdx of the second
%   dimension of res.nmse instead (e.g., a training size).

if nargin < 3 || isempty(sIdx)
    sIdx = find(abs(res.snrDB - snrDB) < 1e-9, 1);
end
for n = 1:size(res.nmse, 4)
    fprintf('\n  %s (column %d)\n  %-20s %16s %16s %6s\n', res.netLabels{n}, sIdx, ...
            'method', 'NMSE [dB]', 'worst [dB]', 'id.');
    for e = 1:numel(res.keys)
        v = squeeze(res.nmse(e, sIdx, :, n));
        if all(isnan(v)), continue; end
        w = squeeze(res.worst(e, sIdx, :, n));
        v = v(~isnan(v)); w = w(~isnan(w));
        [m1, l1, h1] = ut_mc_summary(v(:)', 2);
        [m2, l2, h2] = ut_mc_summary(w(:)', 2);
        idr = squeeze(res.id(e, sIdx, :, n));
        fprintf('  %-20s %8.2f (+-%4.2f) %8.2f (+-%4.2f) %6s\n', res.keys{e}, m1, (h1 - l1) / 2, ...
                m2, (h2 - l2) / 2, num2str(mean(idr(~isnan(idr))), '%.2f'));
    end
end
end
