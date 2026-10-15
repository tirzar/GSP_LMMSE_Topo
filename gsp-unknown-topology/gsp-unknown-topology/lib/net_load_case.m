function mpc = net_load_case(caseName)
%NET_LOAD_CASE  Load a MATPOWER-format test case from the bundled 'cases' folder.
%
%   mpc = NET_LOAD_CASE(caseName) evaluates cases/<caseName>.m. If the file is
%   not bundled but MATPOWER is on the path, MATPOWER's LOADCASE is used.
%
%   Bundled: case14, case24_ieee_rts, case_ieee30 (IEEE 30-bus; MATPOWER's
%   'case30' is a different system), case39, case57, case118.

caseDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'cases');
caseFile = fullfile(caseDir, [caseName '.m']);
if exist(caseFile, 'file') == 2
    oldPath = addpath(caseDir);
    cleaner = onCleanup(@() path(oldPath)); %#ok<NASGU>
    mpc = feval(caseName);
elseif exist('loadcase', 'file') == 2
    mpc = loadcase(caseName);
else
    error('net_load_case:notFound', 'Case "%s" not found in %s.', caseName, caseDir);
end
end
