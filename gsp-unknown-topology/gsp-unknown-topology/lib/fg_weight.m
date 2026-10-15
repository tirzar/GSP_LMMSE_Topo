function [lw, ms] = fg_weight(st, cfg)
%FG_WEIGHT  Line width and marker size of a method style (emphasized or not).
if st.emph
    lw = 1.6 * cfg.fig.lineWidth; ms = cfg.fig.markerSize + 0.5;
else
    lw = 0.85 * cfg.fig.lineWidth; ms = cfg.fig.markerSize;
end
end
