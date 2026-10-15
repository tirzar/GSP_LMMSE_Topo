function cfg = sim_config(mode)
%SIM_CONFIG  Parameters of the simulations of Section VI.
%
%   cfg = SIM_CONFIG() or SIM_CONFIG('full') returns the settings of the paper.
%   cfg = SIM_CONFIG('quick') returns reduced sizes for a short test run (the
%   resulting numbers are not those of the paper).
%
%   Conventions
%     SNR [dB] = -10 log10(sigma^2), with sigma in p.u. on the 100-MVA base,
%       common to all networks (sigma = 0.01 p.u. = 1 MW).
%     NMSE = sum_m p_m MSE_m / sum_m p_m E_m||x - mu_x||^2 (mixture mean mu_x);
%       worst-topology NMSE = max_m MSE_m / E_m||x - mu_x||^2.
%     Each Monte Carlo realization redraws the candidate family, the training
%       data, and the test data.

if nargin < 1 || isempty(mode), mode = 'full'; end
cfg.mode = lower(mode);
cfg.root = fileparts(mfilename('fullpath'));
cfg.dirs.lib     = fullfile(cfg.root, 'lib');
cfg.dirs.cases   = fullfile(cfg.root, 'cases');
cfg.dirs.results = fullfile(cfg.root, 'results');
cfg.dirs.figures = fullfile(cfg.root, 'figures');
cfg.dirs.tables  = fullfile(cfg.root, 'tables');
cfg.isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
cfg.seed = 2026;

% ---------------- networks, candidate families, data models ----------------
cfg.networks  = {'case118', 'case39'};                 % Sections VI-B to VI-E
cfg.netLabels = {'IEEE 118-bus', 'New England 39-bus'};
cfg.M = 20;                                            % number of candidate topologies
cfg.family = struct('pRemove', 0.01, 'nAddMean', 1);   % see NET_TOPOLOGY_FAMILY
cfg.beta = 0.1;                                        % M1: x | m ~ N(0, beta L_m^+)
cfg.pf = struct('sigLoad', 0.1, 'kappa', [0.8 1.2]);   % M2: AC power-flow states
cfg.snrDB    = 0:10:40;                                % SNR grid (M1)
cfg.snr0DB   = 20;                                     % default SNR, M1
cfg.snrM2DB  = 40;                                     % default SNR, M2
cfg.T        = 500;                                    % training samples per topology
cfg.Ttest    = 500;                                    % test samples per topology

% ---------------- DNN benchmark ----------------
cfg.dnn = struct('hidden', [256 256], 'lr', 1e-3, 'batch', 128, 'wd', 1e-5, ...
                 'valFrac', 0.15, 'patience', 15, 'maxEpochs', 200);

% ---------------- Section VI-B: price of an unknown topology (Fig. 1, Table V) ----------------
cfg.B.nMC  = 20;
cfg.B.keys = {'aware_gsp', 'aware_lmmse_ridge', 'blind_lmmse', 'blind_lmmse_ridge', ...
              'mb', 'r1', 'sr', 'mm_gsp_hard', 'mm_gsp_soft', ...
              'nominal_gsp', 'expected_gsp', 'minimax_mse', 'minimax_regret', ...
              'mm_full_hard', 'mm_full_soft', 'dnn'};
cfg.B.figKeys = {'aware_gsp', 'blind_lmmse', 'mb', 'mm_gsp_soft', 'mm_gsp_hard', ...
                 'nominal_gsp', 'expected_gsp', 'minimax_regret', 'mm_full_soft', 'dnn'};

% ---------------- Section VI-C: number of bases and perturbation size (Figs. 2, 3) ----------------
cfg.C.nMC     = 20;
cfg.C.K       = [1 2 5 10 20];     % number of candidate bases in the multi-basis estimator
cfg.C.nSub    = 20;                % random subsets of bases per K (all singletons for K = 1)
cfg.C.scales  = [0.25 0.5 1 2 4];  % perturbation size: (pRemove, nAddMean) x scale
cfg.C.keys    = {'aware_gsp', 'blind_lmmse', 'mb', 'r1', 'sr', 'sb_nom', ...
                 'nominal_gsp', 'expected_gsp', 'minimax_regret', 'mm_gsp_soft', 'mm_gsp_hard'};

% ---------------- Section VI-D: finite training data (Fig. 4) ----------------
cfg.D.nMC  = 20;
cfg.D.T    = [20 50 100 200 500 1000 2000];
cfg.D.keys = {'aware_gsp', 'aware_lmmse_ridge', 'blind_lmmse_ridge', 'mb', 'r1', ...
              'sr', 'expected_gsp', 'mm_gsp_soft', 'mm_full_soft', 'dnn', 'mm_gsp_hard'};
cfg.D.dnnT = [50 200 500 2000];    % training sizes at which the DNN is trained

% ---------------- Section VI-E: P snapshots sharing the topology (Fig. 5) ----------------
cfg.E.nMC   = 20;
cfg.E.T     = [50 500];            % training samples per topology
cfg.E.P     = [1 2 5 10 20 50];    % number of snapshots
cfg.E.Ttest = 1000;
cfg.E.keys  = {'aware_gsp', 'mb', 'mm_gsp_hard', 'mm_gsp_soft', 'mm_full_hard', 'mm_full_soft'};

% ---------------- Section VI-F: transfer across network sizes (Fig. 6, Table VI) ----------------
cfg.F.networks  = {'case14', 'case24_ieee_rts', 'case_ieee30', 'case39', 'case57', 'case118'};
cfg.F.netLabels = {'IEEE 14', 'IEEE RTS 24', 'IEEE 30', 'New England 39', 'IEEE 57', 'IEEE 118'};
cfg.F.nMC     = 5;                 % independent data sets
cfg.F.T       = 20000;             % training samples per network
cfg.F.Ttest   = 5000;              % test samples per network
cfg.F.snrDB   = 20;                % M1; M2 uses cfg.snrM2DB
cfg.F.Q       = 60;                % Chebyshev order
cfg.F.muCheb  = 1e-7;              % stability penalty mu * sum_q q^4 alpha_q^2
cfg.F.poles   = logspace(-2, 2, 4);% rational dictionary 1/(lambda + c_q)
cfg.F.muRat   = 1e-8;              % ridge penalty of the rational dictionary
cfg.F.T0      = [2 5 10 20 50 100 200];   % labeled target samples (few-shot)
cfg.F.nFew    = 50;                % random target training sets per T0

% ---------------- figures ----------------
cfg.fig = struct('visible', 'off', 'fontSize', 8, 'font', 'Times', 'lineWidth', 1, ...
                 'markerSize', 4, 'showCI', true, 'col1', 3.5, 'col2', 7.16);

if strcmp(cfg.mode, 'quick')
    cfg.M = 6;
    cfg.T = 200; cfg.Ttest = 200;
    cfg.snrDB = [10 20 30];
    cfg.dnn.maxEpochs = 5;
    cfg.B.nMC = 2;
    cfg.C.nMC = 2; cfg.C.K = [1 2 6]; cfg.C.nSub = 3; cfg.C.scales = [0.5 1 2];
    cfg.D.nMC = 2; cfg.D.T = [20 100 300]; cfg.D.dnnT = 100;
    cfg.E.nMC = 2; cfg.E.P = [1 5 20]; cfg.E.Ttest = 200; cfg.E.T = [50 200];
    cfg.F.nMC = 2; cfg.F.T = 2000; cfg.F.Ttest = 1000; cfg.F.T0 = [10 50]; cfg.F.nFew = 3;
end

if cfg.isOctave
    cfg.fig.font = 'Liberation Serif';
    if any(strcmp(available_graphics_toolkits(), 'qt')) && ~isempty(getenv('DISPLAY'))
        graphics_toolkit('qt');
    end
    warning('off', 'Octave:singular-matrix');
    warning('off', 'Octave:nearly-singular-matrix');
end
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
addpath(cfg.dirs.lib, cfg.dirs.cases);
for d = {cfg.dirs.results, cfg.dirs.figures, cfg.dirs.tables}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end
end
