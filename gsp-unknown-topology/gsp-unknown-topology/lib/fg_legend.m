function lg = fg_legend(ax, h, labels, cfg, loc, nCol)
%FG_LEGEND  Compact legend (multi-column when supported).
%
%   lg = FG_LEGEND(ax, h, labels, cfg, loc, nCol)

if nargin < 5 || isempty(loc), loc = 'best'; end
if nargin < 6 || isempty(nCol), nCol = 1; end
lg = legend(ax, h, labels, 'Location', loc, 'FontSize', max(cfg.fig.fontSize - 1, 7), 'Box', 'off');
try
    set(lg, 'NumColumns', nCol);
catch
end
try
    set(lg, 'ItemTokenSize', [18 8]);     % MATLAB: shorter line samples
catch
end
end
