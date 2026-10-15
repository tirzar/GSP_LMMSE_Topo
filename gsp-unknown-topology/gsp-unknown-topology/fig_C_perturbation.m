function fig_C_perturbation(res, cfg, mdl)
%FIG_C_PERTURBATION  Fig. 3: NMSE versus the size of the topology perturbation.
%
%   FIG_C_PERTURBATION(res, cfg, mdl) plots the NMSE of the estimators
%   res.keys versus the average number of changed lines per candidate, for
%   the two networks of res (left, right) and data model mdl (default 1, M1).
%   Saves figures/fig_C_perturbation.pdf and .png.

if nargin < 3 || isempty(mdl), mdl = 1; end
cfg = fg_small(cfg);
fig = fg_new(cfg, cfg.fig.col1, 3.3);
pos = {[0.13 0.51 0.355 0.41], [0.625 0.51 0.355 0.41]};
cfgP = cfg; cfgP.fig.markerSize = cfg.fig.markerSize - 0.5;
for n = 1:min(2, numel(res.nets))
    ax = fg_axes(fig, pos{n}, cfg);
    x = mean(res.pt.changes(:, :, n), 2, 'omitnan')';
    h = []; lab = {};
    for e = 1:numel(res.keys)
        V = reshape(res.pt.nmse(e, :, :, n, mdl), numel(x), []);
        [m, lo, hi] = ut_mc_rows(V);
        h(end + 1) = fg_plot(ax, x, m, res.keys{e}, cfgP, lo, hi); %#ok<AGROW>
        lab{end + 1} = ut_short(res.keys{e}); %#ok<AGROW>
    end
    tl = arrayfun(@(v) sprintf('%.1f', v), x, 'UniformOutput', false);
    if numel(tl) > 2, tl(2) = {''}; end            % the two smallest values are close
    set(ax, 'XScale', 'log', 'XTick', x, 'XTickLabel', tl);
    xlim(ax, [min(x) * 0.9, max(x) * 1.1]);
    xlabel(ax, {'average changed lines'; 'per candidate'}); ylabel(ax, 'NMSE [dB]');
    title(ax, sprintf('%s, %s', strrep(strrep(res.netLabels{n}, 'New England 39-bus', 'NE 39'), ...
          'IEEE 118-bus', 'IEEE 118'), res.models{mdl}), 'FontWeight', 'normal');
    if n == 1, h1 = h; lab1 = lab; end
end
lg = fg_legend(ax, h1, lab1, cfg, 'south', 2);
set(lg, 'Position', [0.02 0.025 0.96 0.29]);
fg_export(fig, 'fig_C_perturbation', cfg);
end
