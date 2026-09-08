function run_analyses(name)
%% =========================================================================
%  run_analyses.m  --  MAIN ENTRY POINT 3 of 4
%
%  Dispatcher for every secondary analysis (Monte Carlo, tumbling
%  validation, accommodation studies, numerical error analysis, orbital
%  decay / lifetime). Each analysis keeps its own self-contained physics
%  in lib/analyses/<file>.m -- the small differences between them
%  (rho scaling, wrapTo180 in the SAM binning, ode45 vs ode113, equatorial
%  vs Sun-synchronous, parfor) are deliberate and preserved.
%
%  Usage:
%     run_analyses                       % list the available analyses
%     run_analyses montecarlo_multiparam % run one of them
%
%  Notes:
%    * The working directory is switched to the repo root for the run and
%      restored afterwards, so the scripts' load()/save()/exportgraphics()
%      calls resolve against the repo root as before.
%    * Result variables live in this function's scope only; every analysis
%      also writes its .mat / .png outputs to disk. To keep the workspace
%      variables, open lib/analyses/<file>.m and press Run instead.
%    * montecarlo_barchart, sensitivity_theta0 and sensitivity_isolated
%      load mc_workspace.mat -- run montecarlo_multiparam first.
%    * sensitivity_isolated uses parfor (Parallel Computing Toolbox).
% =========================================================================

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'lib')));
addpath(genpath(fullfile(here, 'adbsat_processed')));

analyses = {
    'montecarlo_multiparam'   'Monte Carlo over theta0, dphi0, rho (phi=75, 450 km); writes mc_workspace.mat'
    'montecarlo_sensitivity'  'Monte Carlo sensitivity, full theta0 range, ode113 (phi=85, 450 km)'
    'montecarlo_barchart'     'Per-geometry Cd_bar with MC error bars (needs mc_workspace.mat)'
    'accom_pdf'               'Residence-time PDFs vs accommodation coefficient (phi=75, all altitudes)'
    'accom_transient'         'Effect of accommodation coefficient on the pitch transient'
    'sensitivity_theta0'      'Isolated initial-pitch-angle sweep (needs mc_workspace.mat)'
    'sensitivity_isolated'    'Isolated theta0 / dphi0 / rho tests, 3 densities (needs mc_workspace.mat; parfor)'
    'tumbling_multiday'       'Multi-day pitch propagation, drag-sail inertia'
    'tumbling_sphere'         'Multi-day pitch propagation, Roberts sphere inertia (validation)'
    'error_convergence'       'n_steps, bin-width, dphi0 and density-scaling convergence checks'
    'decay_map'               'Decay-time contour maps over (alpha, phi) (needs decay_times.mat)'
    'lifetime_kinghele'       'King-Hele orbital lifetime vs altitude, incl. solar-activity band'
};

if nargin == 0 || isempty(name) || any(strcmpi(name, {'list','help','-h','--help'}))
    fprintf('\n  run_analyses <name>\n\n');
    for i = 1:size(analyses,1)
        fprintf('    %-24s  %s\n', analyses{i,1}, analyses{i,2});
    end
    fprintf('\n');
    return
end

name = char(name);
row  = find(strcmpi(name, analyses(:,1)), 1);
if isempty(row)
    error('run_analyses:unknownAnalysis', ...
        'Unknown analysis "%s". Run run_analyses with no arguments for the list.', name);
end

target = fullfile(here, 'lib', 'analyses', [analyses{row,1} '.m']);
old    = cd(here);
restore = onCleanup(@() cd(old));
fprintf('=== run_analyses: %s ===\n', analyses{row,1});
run(target);
end
