function tab_B_comparison(res, cfg)
%TAB_B_COMPARISON  Table V: estimator comparison at the default SNR.
%
%   TAB_B_COMPARISON(res, cfg) writes tables/tab_B_comparison.tex from the
%   results of SIM_B_PRICE: mean and worst-topology NMSE [dB] per network and
%   the design time on the first network. The caption reports the largest
%   95% confidence half-width of the mean NMSE.

s = find(abs(res.snrDB - res.snr0DB) < 1e-9, 1);
nN = numel(res.nets);
groups = {
  'Topology aware (benchmarks)', {'aware_gsp', 'aware_lmmse_ridge'};
  'Topology blind, linear', {'blind_lmmse', 'mb', 'sr', 'r1', 'nominal_gsp', 'expected_gsp', ...
                             'minimax_mse', 'minimax_regret'};
  'Topology blind, nonlinear (data dependent)', {'mm_gsp_hard', 'mm_gsp_soft', 'mm_full_hard', ...
                             'mm_full_soft', 'dnn'}};
rows = {}; hw = 0;
for g = 1:size(groups, 1)
    rows{end + 1} = sprintf('\\multicolumn{%d}{@{}l}{\\emph{%s}}', 2 + 2 * nN, groups{g, 1}); %#ok<AGROW>
    for key = groups{g, 2}
        e = find(strcmp(res.keys, key{1}));
        if isempty(e) || all(isnan(reshape(res.nmse(e, s, :, :), 1, []))), continue; end
        nm = ut_short(key{1});
        info = ut_methods(key{1});
        if strcmp(info.group, 'proposed'), nm = ['\textbf{' nm '}']; end
        line = nm;
        for n = 1:nN
            [m1, l1, h1] = ut_mc_summary(squeeze(res.nmse(e, s, :, n))', 2);
            m2 = ut_mc_summary(squeeze(res.worst(e, s, :, n))', 2);
            hw = max(hw, (h1 - l1) / 2);
            line = sprintf('%s & $%.2f$ & $%.2f$', line, m1, m2);
        end
        rows{end + 1} = sprintf('%s & $%.2g$', line, mean(squeeze(res.tDesign(e, s, :, 1)))); %#ok<AGROW>
    end
end
hdr1 = ' ';
for n = 1:nN
    hdr1 = sprintf('%s& \\multicolumn{2}{c}{%s} ', hdr1, strrep(res.netLabels{n}, 'New England', 'NE'));
end
hdr1 = [hdr1 '& Time'];
hdr2 = ['Estimator' repmat(' & Mean & Worst', 1, nN) ' & [s]'];
cap = sprintf(['Estimator comparison under topology uncertainty (M1, $%g$\\,dB SNR, $M=%d$, uniform prior): ', ...
    'mean and worst-topology NMSE [dB] over $%d$ realizations (95\\%% confidence half-widths of the mean ', ...
    '${\\le}\\,%.2f$\\,dB), and design time on the %s system. Proposed estimators are in bold.'], ...
    res.snr0DB, res.M, size(res.nmse, 3), hw, res.netLabels{1});
tab_write(cfg, 'tab_B_comparison', ['@{}l' repmat('cc', 1, nN) 'c@{}'], {hdr1, hdr2}, rows, cap, ...
          'tab:comparison', false);
end
