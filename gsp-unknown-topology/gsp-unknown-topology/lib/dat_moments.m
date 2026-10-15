function S = dat_moments(fam, D, sigma2, T)
%DAT_MOMENTS  Sample moments used by all estimator designs (Sec. V).
%
%   S = DAT_MOMENTS(fam, D, sigma2, T) uses the first T training samples of each
%   candidate topology, with y = g(L_m, x) + sigma*w, and returns
%     per topology m (conditional moments, (9) and (65)):
%       S.mux(:,m), S.muy(:,m), S.Cxy{m}, S.Cyy{m}, S.trCxx(m),
%       S.dxy(:,m), S.dyy(:,m)   graph-frequency moments in the basis V_m
%     mixture (prior-weighted, (58)-(61)):
%       S.muxMix, S.muyMix, S.CxyMix, S.CyyMix, S.trCxxMix
%     S.Xraw{m}, S.Yraw{m}      the training data (ridge designs, DNN)

M = fam.M;
N = fam.N;
p = fam.p(:);
S.p = p;
S.T = T;
S.sigma2 = sigma2;
S.mux = zeros(N, M); S.muy = zeros(N, M);
S.Cxy = cell(1, M); S.Cyy = cell(1, M);
S.trCxx = zeros(1, M);
S.dxy = zeros(N, M); S.dyy = zeros(N, M);
S.Xraw = cell(1, M); S.Yraw = cell(1, M);
for m = 1:M
    X = D.Xtr{m}(:, 1:T);
    Y = D.Gtr{m}(:, 1:T) + sqrt(sigma2) * D.Wtr{m}(:, 1:T);
    S.Xraw{m} = X; S.Yraw{m} = Y;
    S.mux(:, m) = mean(X, 2);
    S.muy(:, m) = mean(Y, 2);
    Xc = X - S.mux(:, m);
    Yc = Y - S.muy(:, m);
    S.Cxy{m} = Xc * Yc' / T;
    S.Cyy{m} = Yc * Yc' / T;
    S.trCxx(m) = sum(Xc(:) .^ 2) / T;
    Vm = fam.grd{m}.V;
    S.dxy(:, m) = sum(Vm .* (S.Cxy{m} * Vm), 1)';
    S.dyy(:, m) = sum(Vm .* (S.Cyy{m} * Vm), 1)';
end
S.muxMix = S.mux * p;
S.muyMix = S.muy * p;
S.CxyMix = zeros(N); S.CyyMix = zeros(N);
S.trCxxMix = 0;
for m = 1:M
    dx = S.mux(:, m) - S.muxMix;
    dy = S.muy(:, m) - S.muyMix;
    S.CxyMix = S.CxyMix + p(m) * (S.Cxy{m} + dx * dy');
    S.CyyMix = S.CyyMix + p(m) * (S.Cyy{m} + dy * dy');
    S.trCxxMix = S.trCxxMix + p(m) * (S.trCxx(m) + dx' * dx);
end
S.CyyMix = (S.CyyMix + S.CyyMix') / 2;
end
