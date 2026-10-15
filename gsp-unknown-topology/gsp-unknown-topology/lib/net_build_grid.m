function grd = net_build_grid(mpc)
%NET_BUILD_GRID  Admittance matrix, graph Laplacian, and graph Fourier basis.
%
%   grd = NET_BUILD_GRID(mpc) returns
%     .N, .nLines      number of buses and of in-service branches
%     .Ybus            N x N complex bus admittance matrix (pi model with taps,
%                      phase shifters, line charging, and bus shunts; identical to
%                      MATPOWER's makeYbus)
%     .W               weighted adjacency, w_ij = Im{Y_ij} (series susceptance)
%     .L               Laplacian L = diag(W 1) - W; with these weights, L is the
%                      Jacobian of the active-power injections w.r.t. the phase
%                      angles at a flat voltage profile
%     .V, .lambda      eigenvectors and ascending eigenvalues of L (lambda(1) = 0)
%     .ref             index of the reference (slack) bus

BUS_I = 1; BUS_TYPE = 2; GS = 5; BS = 6;
F_BUS = 1; T_BUS = 2; BR_R = 3; BR_X = 4; BR_B = 5; TAP = 9; SHIFT = 10; BR_STATUS = 11;

bus = mpc.bus;
N = size(bus, 1);
e2i = zeros(max(bus(:, BUS_I)), 1);
e2i(bus(:, BUS_I)) = (1:N)';
br = mpc.branch(mpc.branch(:, BR_STATUS) > 0, :);
f = e2i(br(:, F_BUS));
t = e2i(br(:, T_BUS));
nl = size(br, 1);

Ys = 1 ./ (br(:, BR_R) + 1j * br(:, BR_X));
tapRatio = ones(nl, 1);
k = br(:, TAP) ~= 0;
tapRatio(k) = br(k, TAP);
tau = tapRatio .* exp(1j * pi / 180 * br(:, SHIFT));
Ytt = Ys + 1j * br(:, BR_B) / 2;
Yff = Ytt ./ (tau .* conj(tau));
Yft = -Ys ./ conj(tau);
Ytf = -Ys ./ tau;
Ysh = (bus(:, GS) + 1j * bus(:, BS)) / mpc.baseMVA;
Ybus = sparse([f; f; t; t], [f; t; f; t], [Yff; Yft; Ytf; Ytt], N, N) + sparse(1:N, 1:N, Ysh, N, N);

W = full(imag(Ybus));
W(1:N+1:end) = 0;
W = (W + W') / 2;
L = diag(sum(W, 2)) - W;
L = (L + L') / 2;
[V, D] = eig(L);
[lambda, ord] = sort(real(diag(D)), 'ascend');
V = V(:, ord);
lambda(1) = 0;
V(:, 1) = ones(N, 1) / sqrt(N);                 % exact DC eigenvector

ref = find(bus(:, BUS_TYPE) == 3, 1);
if isempty(ref), ref = 1; end
grd = struct('N', N, 'nLines', nl, 'Ybus', Ybus, 'W', W, 'L', L, 'V', V, ...
             'lambda', lambda, 'ref', ref, 'baseMVA', mpc.baseMVA);
end
