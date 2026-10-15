function fig_B_price(res, cfg)
%FIG_B_PRICE  Fig. 1: NMSE versus SNR under topology uncertainty.
%
%   FIG_B_PRICE(res, cfg) plots the Monte Carlo mean NMSE of the estimators
%   cfg.B.figKeys for the two networks of res (SIM_B_PRICE), in two panels
%   stacked vertically with a shared legend, and 95% confidence intervals
%   when cfg.fig.showCI is true. Saves figures/fig_B_price.pdf and .png.

cfg = fg_small(cfg);
cfg.fig.fontSize = 10;
LW = 1.5; MS = 5;
keys = cfg.B.figKeys;
fig = fg_new(cfg, cfg.fig.col1, 6.4);
pos = {[0.17 0.665 0.80 0.29], [0.17 0.285 0.80 0.29]};
hAll = []; lab = {};
for n = 1:numel(res.nets)
    ax = fg_axes(fig, pos{n}, cfg);
    for e = 1:numel(keys)
        k = find(strcmp(res.keys, keys{e}));
        if isempty(k), continue; end
        V = reshape(res.nmse(k, :, :, n), numel(res.snrDB), []);
        [m, lo, hi] = ut_mc_rows(V);
        h = fg_plot(ax, res.snrDB, m, keys{e}, cfg, lo, hi);
        set(h, 'LineWidth', LW, 'MarkerSize', MS);
        if n == 1
            hAll(end + 1) = h; %#ok<AGROW>
            lab{end + 1} = ut_short(keys{e}); %#ok<AGROW>
        end
    end
    xlim(ax, [min(res.snrDB) max(res.snrDB)]);
    set(ax, 'XTick', res.snrDB);
    if n == numel(res.nets), xlabel(ax, 'SNR [dB]'); end
    ylabel(ax, 'NMSE [dB]');
    title(ax, res.netLabels{n}, 'FontWeight', 'normal');
end
legOrder = {'aware_gsp', 'mm_gsp_soft', 'mm_gsp_hard', 'nominal_gsp', 'mm_full_soft', ...
            'expected_gsp', 'mb', 'blind_lmmse', 'minimax_regret', 'dnn'};
shown = keys(ismember(keys, res.keys));
[~, ord] = ismember([legOrder, setdiff(shown, legOrder, 'stable')], shown);
ord = ord(ord > 0);
lg = fg_legend(ax, hAll(ord), lab(ord), cfg, 'south', 2);
set(lg, 'FontSize', cfg.fig.fontSize, 'Position', [0.0 0.0 1.0 0.20]);
fg_export(fig, 'fig_B_price', cfg);
end
