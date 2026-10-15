function fg_export(fig, name, cfg)
%FG_EXPORT  Save a figure as vector PDF and 300-dpi PNG in cfg.dirs.figures.

base = fullfile(cfg.dirs.figures, name);
if exist('exportgraphics', 'file') == 2           % MATLAB R2020a and later
    exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
    exportgraphics(fig, [base '.png'], 'Resolution', 300);
else
    print(fig, [base '.pdf'], '-dpdf');
    print(fig, [base '.png'], '-dpng', '-r300');
end
if ~strcmp(cfg.fig.visible, 'on'), close(fig); end
fprintf('  figure saved: %s.pdf\n', base);
end
