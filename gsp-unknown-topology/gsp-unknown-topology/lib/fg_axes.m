function ax = fg_axes(fig, pos, cfg)
%FG_AXES  Axes with the common style (box, grid, ticks) at normalized pos.

ax = axes('Parent', fig, 'Position', pos);
set(ax, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'FontSize', cfg.fig.fontSize, ...
    'FontName', cfg.fig.font, 'LineWidth', 0.5, 'TickDir', 'in');
hold(ax, 'on');
end
