function tab_A_networks(res, cfg)
%TAB_A_NETWORKS  Table III: test systems.
%
%   TAB_A_NETWORKS(res, cfg) writes tables/tab_A_networks.tex from the
%   results of SIM_A_SETUP: N, number of lines, lambda_2, lambda_max, and the
%   RMS phase angle under data models M1 and M2.

J = numel(res.nets);
rows = cell(1, J);
for j = 1:J
    rows{j} = sprintf('%s & %d & %d & %.2f & %s & %.1f & %.1f', res.netLabels{j}, res.N(j), ...
                      res.nLines(j), res.lam2(j), thousands(res.lamMax(j)), res.rmsM1(j), res.rmsM2(j));
end
cap = ['Test systems: buses, transmission lines, extreme nonzero graph frequencies, ', ...
       'and the root-mean-square (RMS) phase (state) angle, $\sqrt{\mathrm{E}[\|\mathbf{x}\|^2]/N}$, ', ...
       'under models M1 and M2.'];
tab_write(cfg, 'tab_A_networks', '@{}lrrrrrr@{}', ...
    {'& & & & & \multicolumn{2}{c}{RMS angle [$^\circ$]}', ...
     'System & $N$ & Lines & $\lambda_2$ & $\lambda_{\max}$ & M1 & M2'}, ...
    rows, cap, 'tab:networks', false);
end

function s = thousands(v)
s = sprintf('%.0f', v);
if numel(s) > 3, s = [s(1:end-3) '{,}' s(end-2:end)]; end
end
