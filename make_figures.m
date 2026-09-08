function make_figures(name)
%% =========================================================================
%  make_figures.m
%
%  Dispatcher for the figure-generation scripts in lib/figures/. Each one
%  loads sam_results.mat (from SAM.m) and writes .fig + .png into the repo
%  root.
%
%  Usage:
%     make_figures                % list the available figure sets
%     make_figures beta_cd_suite  % run one of them
%
%  The working directory is switched to the repo root for the run and
%  restored afterwards.
%
%  NOTE: the 'sam_sweep' and 'sam_extra' sets load attitude_SAM_sweep_FINAL.mat
%  / longrun_SAM_sweep_FINAL.mat, which come from the 15-orbit runs
%  (lib/variants/NEWMASTER.m and a long-run script not included here). They
%  will error without those files.
% =========================================================================

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'lib')));
addpath(genpath(fullfile(here, 'adbsat_processed')));

figs = {
    'beta_cd_suite'      'figs_beta_cd_suite'      'Cd_bar / beta suite: vs phi, accom, altitude; heatmaps, 3-D surfaces, contours, sensitivity (former plotting.m)'
    'sam_sweep'          'figs_sam_sweep'          'Full SAM-sweep plot set incl. 15-orbit SAM vs random vs no-tumble (former attitude_SAM_plots_all_v2.m)'
    'sam_extra'          'figs_sam_extra'          'SAM-sweep plots plus 7 additional figures (former newfigs.m)'
    'adm'               'figs_adm'                 'Aerodynamic database: Cd and Cm vs pitch angle, Cd at zero incidence vs phi (former admplot.m)'
    'pdf'               'figs_pdf'                 'Residence-time PDF scatter plots, linear and log (former pdfplots.m)'
    'sam_vs_uniform'    'figs_sam_vs_uniform'      'SAM-weighted vs uniform-tumbling Cd_bar comparison (former CDCOMPARE.m)'
    'cd_feedback'       'figs_cd_feedback'         'Cd vs pitch angle, raw ADBSat vs reflected/feedback treatment (former fefdbackplot.m)'
    'angular_velocity'  'figs_angular_velocity'    'Pitch-rate time histories, per geometry and zoomed (former angularvel.m)'
    'density_vs_time'   'figs_density_vs_time'     'Orbit-varying atmospheric density vs time (former newplothelin.m)'
};

if nargin == 0 || isempty(name) || any(strcmpi(name, {'list','help','-h','--help'}))
    fprintf('\n  make_figures <name>\n\n');
    for i = 1:size(figs,1)
        fprintf('    %-18s  %s\n', figs{i,1}, figs{i,3});
    end
    fprintf('\n');
    return
end

name = char(name);
row  = find(strcmpi(name, figs(:,1)), 1);
if isempty(row)
    error('make_figures:unknownFigureSet', ...
        'Unknown figure set "%s". Run make_figures with no arguments for the list.', name);
end

target = fullfile(here, 'lib', 'figures', [figs{row,2} '.m']);
old    = cd(here);
restore = onCleanup(@() cd(old));
fprintf('=== make_figures: %s ===\n', figs{row,1});
run(target);
end
