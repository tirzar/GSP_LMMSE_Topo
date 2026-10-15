function info = ut_methods(key)
%UT_METHODS  Display name, group, and plot style of an estimator.
%
%   info = UT_METHODS(key) returns a struct with fields
%     key, name   estimator key and full name
%     group       'benchmark' | 'proposed' | 'competing' | 'ablation'
%     color, line, marker, emph, fill   plot style: proposed estimators are
%                 solid, competing ones dashed or dotted; emphasized curves
%                 (emph) are drawn thicker; fill gives filled markers.
%   info = UT_METHODS() returns the list of keys.

k = [0 0 0]; gry = [0.45 0.45 0.45]; blu = [0 114 178] / 255; grn = [0 158 115] / 255;
sky = [86 180 233] / 255; mag = [0.78 0.08 0.68]; ver = [213 94 0] / 255;
ora = [230 159 0] / 255; brn = [140 86 75] / 255; vio = [117 112 179] / 255;
tea = [0 110 110] / 255;
%      key                  name                                  group        color line  marker emph   fill
T = {'aware_gsp',         'GSP-LMMSE, known topology',          'benchmark', k,   '-',  'none', true,  false
     'aware_lmmse_ridge', 'LMMSE, known topology',              'benchmark', k,   ':',  '+',    false, false
     'blind_lmmse',       'Topology-blind LMMSE',               'benchmark', gry, '-.', 'x',    false, false
     'blind_lmmse_ridge', 'Topology-blind LMMSE, ridge',        'benchmark', gry, '-.', 'x',    false, false
     'mb',                'Multi-basis',                        'proposed',  blu, '-',  'o',    true,  false
     'r1',                'Rank-one',                           'proposed',  grn, '-',  '^',    false, false
     'sr',                'Shared response',                    'proposed',  sky, '-',  's',    false, false
     'sb_nom',            'Single basis V_1',                   'ablation',  sky, ':',  'd',    false, false
     'mm_gsp_hard',       'Selective GSP',                      'proposed',  mag, '-',  '*',    false, false
     'mm_gsp_soft',       'Multiple-model GSP',                 'proposed',  mag, '-',  'p',    true,  false
     'nominal_gsp',       'GSP-LMMSE, nominal topology',        'competing', ver, '--', 'v',    false, false
     'expected_gsp',      'Expected-graph filter',              'competing', ora, '--', '>',    false, false
     'minimax_mse',       'Minimax MSE',                        'competing', brn, '-.', '<',    false, false
     'minimax_regret',    'Minimax regret',                     'competing', brn, '--', 'h',    false, false
     'mm_full_hard',      'Multiple-model, full LMMSE (hard)',  'competing', tea, ':',  'o',    false, true
     'mm_full_soft',      'Multiple-model, full LMMSE',         'competing', tea, '--', 'o',    false, true
     'dnn',               'DNN',                                'competing', vio, ':',  'd',    false, false};
if nargin == 0
    info = T(:, 1)';
    return;
end
i = find(strcmp(T(:, 1), key));
if isempty(i), error('ut_methods:key', 'Unknown estimator key "%s".', key); end
info = struct('key', key, 'name', T{i, 2}, 'group', T{i, 3}, 'color', T{i, 4}, 'line', T{i, 5}, ...
              'marker', T{i, 6}, 'emph', T{i, 7}, 'fill', T{i, 8});
end
