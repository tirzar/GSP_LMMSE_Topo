function fig_D_training(res, cfg, net)
%FIG_D_TRAINING  Fig. 4: NMSE versus the number of training samples per topology.
%
%   FIG_D_TRAINING(res, cfg, net) plots the estimators res.keys for network
%   net (default 1, IEEE 118-bus): data model M1 (left) and M2 (right). The M2
%   panel is zoomed on the topology-blind estimators; the NMSE of the
%   topology-aware benchmarks, which lie below its range, is printed.
%   Saves figures/fig_D_training.pdf and .png (fig_D_training_net<net> for
%   net > 1).

if nargin < 3 || isempty(net), net = 1; end
keys = res.keys;                                 % selective GSP after multiple-model GSP
iS = find(strcmp(keys, 'mm_gsp_soft')); iH = find(strcmp(keys, 'mm_gsp_hard'));
if ~isempty(iS) && ~isempty(iH)
    keys(iH) = [];
    keys = [keys(1:iS), {'mm_gsp_hard'}, keys(iS + 1:end)];
end
cfg = fg_small(cfg);
fig = fg_new(cfg, cfg.fig.col1, 3.4);
pos = {[0.13 0.50 0.355 0.42], [0.625 0.50 0.355 0.42]};
hAll = []; lab = {};
for mdl = 1:2
    ax = fg_axes(fig, pos{mdl}, cfg);
    yBlind = [];
    for e = 1:numel(keys)
        k = find(strcmp(res.keys, keys{e}));
        V = reshape(res.nmse(k, :, :, net, mdl), numel(res.T), []);
        [m, lo, hi] = ut_mc_rows(V);
        ok = ~isnan(m);
        if ~any(ok), continue; end
        h = fg_plot(ax, res.T(ok), m(ok), keys{e}, cfg, lo(ok), hi(ok));
        if ~strncmp(keys{e}, 'aware', 5), yBlind = [yBlind, m(ok)]; end %#ok<AGROW>
        if mdl == 1
            hAll(end + 1) = h; %#ok<AGROW>
            lab{end + 1} = ut_short(keys{e}); %#ok<AGROW>
        end
    end
    fg_logx(ax, res.T);
    set(ax, 'XTick', [20 100 500 2000], 'XTickLabel', {'20', '100', '500', '2000'});
    if mdl == 2 && ~isempty(yBlind)
        ylim(ax, [floor(min(yBlind)) - 0.5, ceil(max(yBlind)) + 0.5]);
        for kk = {'aware_gsp', 'aware_lmmse_ridge'}
            k = find(strcmp(res.keys, kk{1}));
            if isempty(k), continue; end
            fprintf('  M2, %s at T_m = %d: %.1f dB (outside the plotted range)\n', ...
                    ut_short(kk{1}), res.T(end), ut_mc_summary(reshape(res.nmse(k, end, :, net, mdl), 1, []), 2));
        end
    end
    xlabel(ax, 'samples per topology T_m');
    ylabel(ax, 'NMSE [dB]');
    title(ax, sprintf('%s, %g dB', res.models{mdl}, res.snrDB(mdl)), 'FontWeight', 'normal');
end
lg = fg_legend(ax, hAll, lab, cfg, 'south', 2);
set(lg, 'Position', [0.02 0.0 0.96 0.33]);
name = 'fig_D_training';
if net > 1, name = sprintf('%s_net%d', name, net); end
fg_export(fig, name, cfg);
end
