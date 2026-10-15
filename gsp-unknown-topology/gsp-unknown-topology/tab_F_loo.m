function tab_F_loo(res, cfg)
%TAB_F_LOO  Table VI: leave-one-network-out transfer.
%
%   TAB_F_LOO(res, cfg) writes tables/tab_F_loo.tex from the results of
%   SIM_F_TRANSFER: for each held-out network, the NMSE [dB] of its own
%   GSP-LMMSE estimator and the excess MSE [dB] of the spectral functions
%   learned on the other five networks, under data models M1 and M2 (means
%   over the data sets, measured on test data with the EVD-free filters).

J = numel(res.nets);
dM1 = {'cheb_abs', 'cheb_lmax', 'cheb_graphon', 'rational', 'wiener'};
dM2 = {'cheb_abs', 'cheb_lmax', 'cheb_graphon', 'rational'};
rows = cell(1, J);
for j = 1:J
    line = res.netLabels{j};
    for mdl = 1:2
        own = mean(res.loo.own(j, :, mdl)); den = mean(res.loo.den(j, :, mdl));
        line = sprintf('%s & %.2f', line, 10 * log10(own / den));
        if mdl == 1, ds = dM1; else, ds = dM2; end
        for d = 1:numel(ds)
            k = strcmp(res.dictNames, ds{d});
            v = 10 * log10(mean(res.loo.meas(j, k, :, mdl)) / own);
            if abs(v) < 0.005, v = 0; end
            line = sprintf('%s & $%.2f$', line, v);
        end
    end
    rows{j} = line;
end
hdr1 = ['& \multicolumn{6}{c}{M1, SNR $=' sprintf('%g', res.snrDB(1)) '$\,dB} & ' ...
        '\multicolumn{5}{c}{M2, SNR $=' sprintf('%g', res.snrDB(2)) '$\,dB}'];
hdr2 = ['Target network & Own & Cheb. & Cheb.\ $\lambda/\lambda_{\max}$ & Cheb.\ $\lambda/N$ & ' ...
        'Rational & Wiener & Own & Cheb. & Cheb.\ $\lambda/\lambda_{\max}$ & Cheb.\ $\lambda/N$ & Rational'];
hdr3 = ['& NMSE & abs. & \cite{Defferrard2016} & \cite{Ruiz_Chamon_Ribeiro_GSP_2021} & & (lin.) & ' ...
        'NMSE & abs. & \cite{Defferrard2016} & \cite{Ruiz_Chamon_Ribeiro_GSP_2021} & '];
cap = sprintf(['Leave-one-network-out spectral transfer under M1 (SNR $=%g$\\,dB) and M2 (SNR $=%g$\\,dB). ', ...
    'The ``Own'''' columns give the topology-aware GSP-LMMSE NMSE; all other columns give the excess MSE ', ...
    'relative to this benchmark (dB). Chebyshev filters use $Q=%d$. The Wiener baseline uses the fixed ', ...
    'response $\\beta/(\\beta\\lambda+\\sigma^2)$. Results are averaged over %d realizations.'], ...
    res.snrDB(1), res.snrDB(2), res.Q, size(res.loo.own, 2));
tab_write(cfg, 'tab_F_loo', 'l|rrrrrr|rrrrr', {hdr1, hdr2, hdr3}, rows, cap, 'tab:loo', true);
end
