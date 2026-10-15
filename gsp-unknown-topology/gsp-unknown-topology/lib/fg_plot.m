function h = fg_plot(ax, x, yDB, key, cfg, loDB, hiDB)
%FG_PLOT  Plot one method's curve in its registered style (UT_METHODS).
%
%   h = FG_PLOT(ax, x, yDB, key, cfg) plots yDB versus x; with loDB and hiDB,
%   95% confidence intervals are drawn by FG_CI (when cfg.fig.showCI). The
%   emphasized methods (UT_METHODS: known-topology GSP-LMMSE, multi-basis,
%   multiple-model GSP) are drawn thicker, the other ones slightly thinner.

st = ut_methods(key);
[lw, ms] = fg_weight(st, cfg);
h = plot(ax, x, yDB, 'Color', st.color, 'LineStyle', st.line, 'Marker', st.marker, ...
         'LineWidth', lw, 'MarkerSize', ms);
if st.fill, set(h, 'MarkerFaceColor', st.color); end
if nargin >= 7
    fg_ci(ax, x, loDB, hiDB, st.color, cfg);
end
end
