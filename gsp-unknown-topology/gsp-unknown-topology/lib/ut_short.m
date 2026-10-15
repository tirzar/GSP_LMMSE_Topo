function s = ut_short(key)
%UT_SHORT  Short estimator label for legends and tables (full names: UT_METHODS).
switch key
    case 'aware_gsp',          s = 'GSP-LMMSE (known top.)';
    case 'aware_lmmse_ridge',  s = 'LMMSE (known top.)';
    case 'blind_lmmse',        s = 'Blind LMMSE';
    case 'blind_lmmse_ridge',  s = 'Blind LMMSE (ridge)';
    case 'nominal_gsp',        s = 'Nominal-topology GSP';
    case 'expected_gsp',       s = 'Expected graph';
    case 'mm_full_hard',       s = 'Multiple-model, full (hard)';
    case 'mm_full_soft',       s = 'Multiple-model, full';
    otherwise
        info = ut_methods(key);
        s = info.name;
end
end
