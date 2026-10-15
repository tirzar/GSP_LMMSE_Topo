function [ests, tDesign] = est_design_set(keys, fam, S, cfg)
%EST_DESIGN_SET  Design the estimators named in the cell array keys.
%
%   [ests, tDesign] = EST_DESIGN_SET(keys, fam, S, cfg) returns the designed
%   estimators (in the order of keys) and their design times [s], from the
%   candidate family fam (NET_TOPOLOGY_FAMILY) and the sample moments S
%   (DAT_MOMENTS). Keys:
%     proposed    'mb' 'r1' 'sr' 'mm_gsp_hard' 'mm_gsp_soft'
%     ablation    'sb_nom'
%     benchmarks  'aware_gsp' 'aware_lmmse_ridge' 'blind_lmmse' 'blind_lmmse_ridge'
%     competing   'nominal_gsp' 'expected_gsp' 'minimax_mse' 'minimax_regret'
%                 'mm_full_hard' 'mm_full_soft' 'dnn'

ests = cell(1, numel(keys));
tDesign = zeros(1, numel(keys));
for i = 1:numel(keys)
    t0 = tic;
    switch keys{i}
        case 'mb',                ests{i} = est_multibasis(fam, S);
        case 'r1',                ests{i} = est_rank_one(fam, S);
        case 'sr',                ests{i} = est_shared_response(fam, S, fam.p, 'sr');
        case 'sb_nom',            ests{i} = est_single_basis(fam.grd{1}.V, S, 'sb_nom');
        case 'blind_lmmse',       ests{i} = est_blind_lmmse(fam, S, 'none');
        case 'blind_lmmse_ridge', ests{i} = est_blind_lmmse(fam, S, 'ridge');
        case 'aware_gsp',         ests{i} = est_bank(fam, S, 'gsp');
        case 'aware_lmmse_ridge', ests{i} = est_bank(fam, S, 'lmmse_ridge');
        case 'nominal_gsp'                      % GSP-LMMSE of the nominal topology, from D_1 only
            V = fam.grd{1}.V;
            A = V * ((S.dxy(:, 1) ./ S.dyy(:, 1)) .* V');
            ests{i} = est_affine('nominal_gsp', A, S.mux(:, 1), S.muy(:, 1), fam.N);
        case 'expected_gsp'                     % graph filter of sum_m p_m L_m
            Lbar = zeros(fam.N);
            for m = 1:fam.M, Lbar = Lbar + fam.p(m) * fam.grd{m}.L; end
            [Vbar, Dbar] = eig((Lbar + Lbar') / 2);
            [~, o] = sort(diag(Dbar));
            ests{i} = est_single_basis(Vbar(:, o), S, 'expected_gsp');
        case 'minimax_mse',       ests{i} = est_minimax(fam, S, 'mse');
        case 'minimax_regret',    ests{i} = est_minimax(fam, S, 'regret');
        case 'mm_gsp_hard',       ests{i} = est_mm(fam, S, 'gsp', 'hard');
        case 'mm_gsp_soft',       ests{i} = est_mm(fam, S, 'gsp', 'soft');
        case 'mm_full_hard',      ests{i} = est_mm(fam, S, 'full', 'hard');
        case 'mm_full_soft',      ests{i} = est_mm(fam, S, 'full', 'soft');
        case 'dnn',               ests{i} = est_dnn(S, cfg.dnn);
        otherwise
            error('est_design_set:key', 'Unknown estimator key "%s".', keys{i});
    end
    ests{i}.key = keys{i};
    tDesign(i) = toc(t0);
end
end
