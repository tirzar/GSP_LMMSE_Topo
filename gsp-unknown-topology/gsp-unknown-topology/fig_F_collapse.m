function fig_F_collapse(res, cfg)
%FIG_F_COLLAPSE  Fig. 6: GSP-LMMSE responses (55) of the six networks.
%
%   FIG_F_COLLAPSE(res, cfg) plots, versus the graph frequency, the GSP-LMMSE
%   response of every network (mean over the data sets) under data model M1
%   (left), with the linearized-model response beta/(beta lambda + sigma^2),
%   and under M2 (right), with 1/lambda. Saves figures/fig_F_collapse.pdf and .png.

cfg = fg_small(cfg);
J = numel(res.nets);
fig = fg_new(cfg, cfg.fig.col1, 2.35);
cols = lines(J);
mk = {'o', 's', '^', 'd', 'v', '>'};
pos = {[0.13 0.35 0.36 0.55], [0.62 0.35 0.36 0.55]};
for mdl = 1:2
    ax = fg_axes(fig, pos{mdl}, cfg);
    h = []; lab = {};
    lamAll = [];
    for j = 1:J
        A = res.collapse{j, mdl};                % [lambda, f_1, ..., f_R]
        lam = A(2:end, 1); Fr = A(2:end, 2:end); % the DC mode is excluded
        f = mean(Fr, 2);
        lamAll = [lamAll; lam]; %#ok<AGROW>
        h(end + 1) = plot(ax, lam, f, 'LineStyle', 'none', 'Marker', mk{j}, 'Color', cols(j, :), ...
                          'MarkerSize', 2); %#ok<AGROW>
        if size(Fr, 2) > 1
            se = std(Fr, 0, 2) / sqrt(size(Fr, 2));
            fg_ci(ax, lam, max(f - 1.96 * se, 1e-6), f + 1.96 * se, cols(j, :), cfg);
        end
        lab{end + 1} = strrep(strrep(res.netLabels{j}, 'New England', 'NE'), 'IEEE RTS', 'RTS'); %#ok<AGROW>
    end
    lg = logspace(log10(min(lamAll)), log10(max(lamAll)), 200);
    sigma2 = ut_snr(res.snrDB(mdl));
    if mdl == 1
        h(end + 1) = plot(ax, lg, res.beta ./ (res.beta * lg + sigma2), 'k-', 'LineWidth', 0.9);
    else
        h(end + 1) = plot(ax, lg, 1 ./ lg, 'k-', 'LineWidth', 0.9);
    end
    lab{end + 1} = 'reference';
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XTick', [1 10 100 1000], 'YTick', [1e-3 1e-2 1e-1 1 10]);
    xlim(ax, [min(lamAll) / 1.3, max(lamAll) * 1.3]);
    ylim(ax, [5e-4 10]);
    xlabel(ax, 'graph frequency \lambda');
    if mdl == 1, ylabel(ax, 'response f(\lambda)'); end
    title(ax, sprintf('%s, %g dB', res.models{mdl}, res.snrDB(mdl)), 'FontWeight', 'normal');
end
L = fg_legend(ax, h, lab, cfg, 'south', 4);
set(L, 'Position', [0.02 0.005 0.96 0.13]);
fg_export(fig, 'fig_F_collapse', cfg);
end
