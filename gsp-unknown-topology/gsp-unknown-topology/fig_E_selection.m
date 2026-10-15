function fig_E_selection(res, cfg, net)
%FIG_E_SELECTION  Fig. 5: topology selection using P snapshots.
%
%   FIG_E_SELECTION(res, cfg, net) plots, for network net (default 1, IEEE
%   118-bus):
%   left   identification rate of the approximate MAP rule (38) with the
%          graph-stationary surrogate (selective GSP) and with the
%          full-covariance likelihood, for the smallest (dashed) and largest
%          (solid) training sizes in res.T;
%   right  NMSE at the largest training size of the multiple-model estimators
%          (soft weighting and hard selection), with the multi-basis estimator
%          and the topology-aware benchmark.
%   Saves figures/fig_E_selection.pdf and .png (fig_E_selection_net<net> for
%   net > 1).

if nargin < 3 || isempty(net), net = 1; end
cfg = fg_small(cfg);
fig = fg_new(cfg, cfg.fig.col1, 3.45);
pos = {[0.13 0.52 0.355 0.40], [0.625 0.52 0.355 0.40]};
P = res.P;
tStyle = {'--', '-'};
tOrd = [find(res.T == min(res.T), 1), find(res.T == max(res.T), 1)];
lw = cfg.fig.lineWidth; ms = cfg.fig.markerSize;
for panel = 1:2
    ax = fg_axes(fig, pos{panel}, cfg);
    h = []; lab = {};
    if panel == 1
        keys = {'mm_gsp_hard', 'mm_full_hard'};
        for e = 1:numel(keys)
            k = find(strcmp(res.keys, keys{e}));
            st = ut_methods(keys{e});
            for ti = 1:numel(tOrd)
                t = tOrd(ti);
                Iall = reshape(res.id(k, :, t, :, net), numel(P), []);
                I = mean(Iall, 2)';
                h(end + 1) = plot(ax, P, I, 'Color', st.color, 'Marker', st.marker, ...
                    'LineStyle', tStyle{ti}, 'LineWidth', lw, 'MarkerSize', ms); %#ok<AGROW>
                if st.fill, set(h(end), 'MarkerFaceColor', st.color); end
                se = std(Iall, 0, 2)' / sqrt(size(Iall, 2));
                fg_ci(ax, P, max(I - 1.96 * se, 0), min(I + 1.96 * se, 1), st.color, cfg);
                lab{end + 1} = sprintf('%s, T_m = %d', likelihoodName(keys{e}), res.T(t)); %#ok<AGROW>
            end
        end
        h(end + 1) = plot(ax, P, ones(size(P)) / res.M, 'k:', 'LineWidth', 0.75);
        lab{end + 1} = 'chance, 1/M';
        ylim(ax, [0 1]);
        ylabel(ax, 'identification rate');
        ttl = 'Identification';
    else
        keys = {'aware_gsp', 'mm_gsp_soft', 'mm_gsp_hard', 'mm_full_soft', 'mm_full_hard', 'mb'};
        t = tOrd(end);
        for e = 1:numel(keys)
            k = find(strcmp(res.keys, keys{e}));
            if isempty(k), continue; end
            V = reshape(res.nmse(k, :, t, :, net), numel(P), []);
            [m, lo, hi] = ut_mc_rows(V);
            h(end + 1) = fg_plot(ax, P, m, keys{e}, cfg, lo, hi); %#ok<AGROW>
            lab{end + 1} = ut_short(keys{e}); %#ok<AGROW>
        end
        ylabel(ax, 'NMSE [dB]');
        ttl = sprintf('NMSE, T_m = %d', res.T(t));
    end
    fg_logx(ax, P);
    xlabel(ax, 'snapshots P');
    title(ax, ttl, 'FontWeight', 'normal');
    lg = fg_legend(ax, h, lab, cfg, 'south', 2);
    if panel == 1
        set(lg, 'Position', [0.01 0.165 0.98 0.17]);
    else
        set(lg, 'Position', [0.01 0.0 0.98 0.17]);
    end
end
name = 'fig_E_selection';
if net > 1, name = sprintf('%s_net%d', name, net); end
fg_export(fig, name, cfg);
end

function s = likelihoodName(key)
if strncmp(key, 'mm_gsp', 6)
    s = 'Selective GSP';
else
    s = 'Mult.-model, full';
end
end
