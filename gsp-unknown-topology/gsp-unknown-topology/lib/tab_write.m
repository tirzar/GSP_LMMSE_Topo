function tab_write(cfg, name, colSpec, header, rows, caption, label, wide)
%TAB_WRITE  Write a LaTeX table (IEEEtran style) to <cfg.dirs.tables>/<name>.tex.
%
%   TAB_WRITE(cfg, name, colSpec, header, rows, caption, label, wide)
%     colSpec  column specification, e.g., 'lccc'
%     header   cell array of header rows (without the trailing \\)
%     rows     cell array of body rows
%     wide     true for a two-column table* (default false)

if nargin < 8, wide = false; end
env = 'table';
if wide, env = 'table*'; end
fn = fullfile(cfg.dirs.tables, [name '.tex']);
fid = fopen(fn, 'w');
fprintf(fid, '\\begin{%s}[t]\n\\centering\n\\caption{%s}\n\\label{%s}\n', env, caption, label);
fprintf(fid, '\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n\\begin{tabular}{%s}\n\\hline\n', colSpec);
for i = 1:numel(header)
    fprintf(fid, '%s\\\\\n', header{i});
end
fprintf(fid, '\\hline\n');
for i = 1:numel(rows)
    fprintf(fid, '%s\\\\\n', rows{i});
end
fprintf(fid, '\\hline\n\\end{tabular}\n\\end{%s}\n', env);
fclose(fid);
fprintf('  table saved: %s\n', fn);
end
