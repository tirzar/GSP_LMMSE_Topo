function fg_logx(ax, ticks)
%FG_LOGX  Logarithmic x-axis with plain tick labels (e.g., 20, 50, 100).
set(ax, 'XScale', 'log', 'XTick', ticks, 'XTickLabel', arrayfun(@num2str, ticks, 'UniformOutput', false));
xlim(ax, [min(ticks), max(ticks)]);
end
