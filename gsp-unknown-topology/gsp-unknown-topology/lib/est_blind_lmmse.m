function est = est_blind_lmmse(fam, S, reg)
%EST_BLIND_LMMSE  Topology-blind LMMSE estimator (LMAMSE), eq. (12).
%
%   est = EST_BLIND_LMMSE(fam, S, reg) returns x_hat = mux + Cxy Cyy^{-1}(y - muy)
%   with the sample mixture moments S. reg = 'none' uses Cyy^+ (eq. (12) with
%   sample moments); reg = 'ridge' adds a cross-validated diagonal loading
%   (EST_RIDGE) computed from the pooled, prior-weighted training data.

switch reg
    case 'none'
        A = S.CxyMix * pinv(S.CyyMix);
        key = 'blind_lmmse';
    case 'ridge'
        Xc = cell2mat(cellfun(@(X) X - S.muxMix, S.Xraw, 'UniformOutput', false));
        Yc = cell2mat(cellfun(@(Y) Y - S.muyMix, S.Yraw, 'UniformOutput', false));
        w = repelem(S.p(:)' / S.T, 1, S.T);
        A = est_ridge(Xc, Yc, w);
        key = 'blind_lmmse_ridge';
end
est = est_affine(key, A, S.muxMix, S.muyMix, fam.N ^ 2);
end
