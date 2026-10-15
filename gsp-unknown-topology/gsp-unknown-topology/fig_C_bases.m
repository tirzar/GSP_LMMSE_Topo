function fig_C_bases(res, cfg)
%FIG_C_BASES  Fig. 2: multi-basis estimator with K of the M candidate bases.
%
%   FIG_C_BASES(res, cfg) plots the excess NMSE of the multi-basis estimator
%   with K bases over the topology-blind LMMSE estimator (12), averaged over
%   the topologies (left) and for the worst topology (right), versus K. Color
%   and marker indicate the network, the line style the data model (M1 solid,
%   M2 dashed). Values are averaged over the subsets of bases and the
%   realizations (SIM_C_HIERARCHY). Saves figures/fig_C_bases.pdf and .png.

cfg = fg_small(cfg);
fig = fg_new(cfg, cfg.fig.col1, 2.3);
pos = {[0.13 0.40 0.355 0.48], [0.625 0.40 0.355 0.48]};
netCol = {[0 114 178] / 255, [213 94 0] / 255}; netMk = {'o', 's'};
mdlLs = {'-', '--'};
fields = {'nmse', 'worst'}; refF = {'ref', 'refWorst'};
ttl = {'Average', 'Worst topology'};
for p = 1:2
    ax = fg_axes(fig, pos{p}, cfg);
    h = []; lab = {};
    for mdl = 1:2
        for n = 1:numel(res.nets)
            V = res.kb.(fields{p})(:, :, :, n, mdl);          % K x subsets x R
            m = zeros(1, numel(res.K));
            for kk = 1:numel(res.K)
                v = V(kk, :, :); m(kk) = mean(v(~isnan(v)));
            end
            refR = reshape(res.kb.(refF{p})(1, :, n, mdl), 1, []);
            y = 10 * log10(m / mean(refR));
            h(end + 1) = plot(ax, res.K, y, 'Color', netCol{n}, 'LineStyle', mdlLs{mdl}, ...
                'Marker', netMk{n}, 'LineWidth', cfg.fig.lineWidth, 'MarkerSize', cfg.fig.markerSize); %#ok<AGROW>
            lab{end + 1} = sprintf('%s, %s', shortNet(res.netLabels{n}), res.models{mdl}); %#ok<AGROW>
            if cfg.fig.showCI
                lo = NaN(size(y)); hi = lo;
                for kk = 1:numel(res.K)
                    mr = squeeze(meanOmitNaN(V(kk, :, :), 2))';   % per-realization mean over subsets
                    [c, l, u] = ut_mc_summary(mr ./ refR, 2);
                    lo(kk) = y(kk) + (l - c); hi(kk) = y(kk) + (u - c);
                end
                fg_ci(ax, res.K, lo, hi, netCol{n}, cfg);
            end
        end
    end
    plot(ax, [min(res.K) max(res.K)], [0 0], 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    fg_logx(ax, res.K);
    xlabel(ax, 'candidate bases K'); ylabel(ax, 'excess over (12) [dB]');
    title(ax, ttl{p}, 'FontWeight', 'normal');
end
lg = fg_legend(ax, h, lab, cfg, 'south', 2);
set(lg, 'Position', [0.05 0.0 0.90 0.17]);
fg_export(fig, 'fig_C_bases', cfg);
end

function m = meanOmitNaN(v, dim)
ok = ~isnan(v); v(~ok) = 0;
m = sum(v, dim) ./ max(sum(ok, dim), 1);
end

function s = shortNet(s)
s = strrep(strrep(s, 'New England 39-bus', 'NE 39'), 'IEEE 118-bus', 'IEEE 118');
end
