function res = sim_A_setup(cfg)
%SIM_A_SETUP  Section VI-A: test systems and candidate families (Table III).
%
%   res = SIM_A_SETUP(cfg) validates the Laplacian of every test system
%   (NET_CHECK_LAPLACIAN) and reports N, the number of in-service lines,
%   lambda_2, lambda_max, and the RMS phase angle under data models M1 and M2.
%   For the candidate families of Sections VI-B to VI-E, it reports the number
%   of changed lines and the basis alignment (1/N) sum_n (v_{1,n}' v_{m,n})^2.
%   Output: results/sim_A.mat and tables/tab_A_networks.tex.

if nargin < 1, cfg = sim_config(); end
fprintf('\n=== Section VI-A: setup (%s mode) ===\n', cfg.mode);
names = cfg.F.networks; labels = cfg.F.netLabels;
J = numel(names);
res = struct('nets', {names}, 'netLabels', {labels});
res.N = zeros(1, J); res.nLines = res.N; res.lam2 = res.N; res.lamMax = res.N;
res.rmsM1 = res.N; res.rmsM2 = res.N; res.jacErr = res.N;
for j = 1:J
    mpc = net_load_case(names{j});
    g = net_build_grid(mpc);
    [ok, qa] = net_check_laplacian(g, mpc);
    if ~ok, error('sim_A_setup:qa', 'Laplacian check failed for %s.', names{j}); end
    ut_seed(cfg.seed + j);
    X1 = dat_prior_states(g, cfg.beta, 2000);
    X2 = dat_pf_states(mpc, g, 500, cfg.pf);
    res.N(j) = g.N; res.nLines(j) = g.nLines; res.lam2(j) = g.lambda(2); res.lamMax(j) = g.lambda(end);
    res.rmsM1(j) = sqrt(mean(X1(:) .^ 2)) * 180 / pi;
    res.rmsM2(j) = sqrt(mean(X2(:) .^ 2)) * 180 / pi;
    res.jacErr(j) = qa.jacobianErr;
    fprintf('  %-16s N = %3d, lines = %3d, lambda_2 = %6.2f, lambda_max = %7.1f, RMS angle M1 %.1f deg, M2 %.1f deg\n', ...
            labels{j}, g.N, g.nLines, g.lambda(2), g.lambda(end), res.rmsM1(j), res.rmsM2(j));
end
nN = numel(cfg.networks); R = max([cfg.B.nMC, cfg.C.nMC, cfg.D.nMC, cfg.E.nMC]);
res.fam.changes = cell(1, nN); res.fam.align = cell(1, nN);
for n = 1:nN
    mpc = net_load_case(cfg.networks{n});
    ch = []; al = [];
    for r = 1:R
        ut_seed(ut_seed_of(cfg, n, r, 1));
        fam = net_topology_family(mpc, cfg.M, cfg.family);
        ch = [ch, fam.nChanges(2:end)]; al = [al, fam.align(2:end)]; %#ok<AGROW>
    end
    res.fam.changes{n} = ch; res.fam.align{n} = al;
    fprintf('  %s, %d families: changed lines %d-%d, basis alignment %.2f-%.2f\n', ...
            cfg.netLabels{n}, R, min(ch), max(ch), min(al), max(al));
end
res.R = R; res.M = cfg.M; res.beta = cfg.beta;
save(fullfile(cfg.dirs.results, 'sim_A.mat'), 'res', '-v7');
tab_A_networks(res, cfg);
end
