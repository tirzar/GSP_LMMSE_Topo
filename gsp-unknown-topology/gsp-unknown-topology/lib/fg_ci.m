function fg_ci(ax, x, lo, hi, color, cfg)
%FG_CI  95% confidence intervals as vertical bars with short caps.
%
%   FG_CI(ax, x, lo, hi, color, cfg) draws, at each x(i), a vertical bar from
%   lo(i) to hi(i) (same units as the plotted curve) with caps whose width is a
%   fixed fraction of the axis (multiplicative on a logarithmic x-axis). Drawn
%   only when cfg.fig.showCI is true; NaN entries are skipped. The bars are
%   excluded from legends.

if ~cfg.fig.showCI, return; end
x = x(:)'; lo = lo(:)'; hi = hi(:)';
isLog = strcmp(get(ax, 'XScale'), 'log') || (all(x > 0) && max(x) / min(x) >= 10);  % axis set to log later
hold(ax, 'on');
color = color + 0.35 * (1 - color);    % lighter than the curve, so that the bars stay secondary
for i = 1:numel(x)
    if any(isnan([x(i) lo(i) hi(i)])) || hi(i) <= lo(i), continue; end
    if isLog
        xl = x(i) / 1.04; xr = x(i) * 1.04;
    else
        w = 0.008 * max(abs(diff(get(ax, 'XLim'))), eps);
        xl = x(i) - w; xr = x(i) + w;
    end
    plot(ax, [x(i) x(i) NaN xl xr NaN xl xr], [lo(i) hi(i) NaN lo(i) lo(i) NaN hi(i) hi(i)], '-', ...
         'Color', color, 'LineWidth', 0.4, 'HandleVisibility', 'off');
end
end
