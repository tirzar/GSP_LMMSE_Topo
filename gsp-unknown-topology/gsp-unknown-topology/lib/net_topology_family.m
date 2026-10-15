function fam = net_topology_family(mpc, M, opts)
%NET_TOPOLOGY_FAMILY  Candidate topologies on a common bus set (Sec. VI-A).
%
%   fam = NET_TOPOLOGY_FAMILY(mpc, M, opts) returns M MATPOWER structs: the
%   nominal network (m = 1) and M-1 variants in which every in-service branch
%   is removed with probability opts.pRemove and Poisson(opts.nAddMean) new
%   lines are added between buses at hop distance two (parameters copied from
%   a random existing line without tap or phase shift). Each variant differs
%   from the nominal network by at least one line, differs from all previous
%   variants, and is accepted only if its Laplacian passes NET_CHECK_LAPLACIAN
%   (connected, distinct eigenvalues, ...); otherwise it is redrawn.
%
%   fam.mpc{m}, fam.grd{m}  case data and graph model of candidate m
%   fam.nChanges(m)         number of removed + added lines w.r.t. nominal
%   fam.align(m)            (1/N) sum_n (v_{1,n}' v_{m,n})^2, basis alignment
%   fam.p                   uniform prior (M x 1)

if ~isfield(opts, 'pRemove'),  opts.pRemove = 0.01; end
if ~isfield(opts, 'nAddMean'), opts.nAddMean = 1; end
if ~isfield(opts, 'maxTries'), opts.maxTries = 10000; end

N = size(mpc.bus, 1);
e2i = zeros(max(mpc.bus(:, 1)), 1);
e2i(mpc.bus(:, 1)) = (1:N)';
onLines = find(mpc.branch(:, 11) > 0);
A = sparse(e2i(mpc.branch(onLines, 1)), e2i(mpc.branch(onLines, 2)), 1, N, N);
A = (A + A') > 0;
A2 = (double(A) * double(A) > 0) & ~A;
A2(1:N+1:end) = false;
[ci, cj] = find(triu(A2));
cand = [mpc.bus(ci, 1), mpc.bus(cj, 1)];
templates = onLines(mpc.branch(onLines, 9) == 0 & mpc.branch(onLines, 10) == 0);

fam.mpc = cell(1, M);
fam.grd = cell(1, M);
fam.nChanges = zeros(1, M);
fam.mpc{1} = mpc;
fam.grd{1} = net_build_grid(mpc);
[ok, ~] = net_check_laplacian(fam.grd{1}, mpc);
if ~ok
    error('net_topology_family:nominal', 'The nominal Laplacian fails the QA.');
end
signatures = {edgeSignature(fam.grd{1})};
m = 2;
tries = 0;
while m <= M
    tries = tries + 1;
    if tries > opts.maxTries
        error('net_topology_family:tries', 'Could not draw %d valid candidates.', M);
    end
    rem = onLines(rand(numel(onLines), 1) < opts.pRemove);
    nAdd = poissrnd_local(opts.nAddMean);
    nAdd = min(nAdd, size(cand, 1));
    if isempty(rem) && nAdd == 0
        continue;
    end
    v = mpc;
    v.branch(rem, 11) = 0;
    add = cand(randperm(size(cand, 1), nAdd), :);
    for k = 1:nAdd
        row = mpc.branch(templates(randi(numel(templates))), :);
        row(1:2) = add(k, :);
        v.branch(end + 1, :) = row; %#ok<AGROW>
    end
    g = net_build_grid(v);
    [ok, ~] = net_check_laplacian(g, v);
    sig = edgeSignature(g);
    if ~ok || any(strcmp(signatures, sig))
        continue;
    end
    fam.mpc{m} = v;
    fam.grd{m} = g;
    fam.nChanges(m) = numel(rem) + nAdd;
    signatures{end + 1} = sig; %#ok<AGROW>
    m = m + 1;
end
fam.M = M;
fam.N = N;
fam.p = ones(M, 1) / M;
fam.align = cellfun(@(g) mean(sum(fam.grd{1}.V .* g.V, 1) .^ 2), fam.grd);
end

function s = edgeSignature(g)
[i, j] = find(triu(g.W) ~= 0);
s = sprintf('%d,', sub2ind(size(g.W), i, j));
end

function k = poissrnd_local(mu)
% Poisson sample by inversion (no Statistics Toolbox needed)
k = 0;
p = exp(-mu);
s = p;
u = rand;
while u > s
    k = k + 1;
    p = p * mu / k;
    s = s + p;
end
end
