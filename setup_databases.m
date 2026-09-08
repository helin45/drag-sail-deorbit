function setup_databases(step, adbsat_results_path)
%% =========================================================================
%  setup_databases.m  --  MAIN ENTRY POINT 1 of 4
%
%  Builds the small input databases the sweep and analyses depend on:
%
%    'inertia'   ->  inertia_tensors.mat   (pitch inertia per apex angle)
%    'geometry'  ->  ref_geometry.mat      (reference area / length per
%                                           apex angle, read from the
%                                           "obj files/" geometries via
%                                           get_ref_geometry)
%    'rename'    ->  copies raw ADBSat run outputs into the
%                    adbsat_processed/<alt>km/<phi>deg_CLL_accom_<a>.mat
%                    layout the rest of the code expects. Requires the
%                    path to the ADBSat 'inou/results' folder as the
%                    second argument.
%    'all'       ->  'inertia' + 'geometry' (not 'rename').
%
%  Usage:
%     setup_databases                      % same as 'all'
%     setup_databases inertia
%     setup_databases geometry
%     setup_databases rename  'D:\ADBSat-master\inou\results'
%
%  inertia_tensors.mat and ref_geometry.mat are already committed to the
%  repo; you only need this to regenerate them or after changing a geometry.
% =========================================================================

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'lib')));
old = cd(here);
restore = onCleanup(@() cd(old));

if nargin == 0 || isempty(step); step = 'all'; end
step = lower(char(step));

switch step
    case 'inertia'
        build_inertia();
    case 'geometry'
        build_geometry(here);
    case 'rename'
        if nargin < 2 || isempty(adbsat_results_path)
            error('setup_databases:missingPath', ...
                ['Provide the ADBSat results folder, e.g.\n' ...
                 '   setup_databases rename ''D:\\ADBSat-master\\inou\\results''']);
        end
        rename_adbsat(here, adbsat_results_path);
    case 'all'
        build_inertia();
        build_geometry(here);
    otherwise
        error('setup_databases:unknownStep', ...
            'Unknown step "%s" (use inertia | geometry | rename | all).', step);
end
end


% -------------------------------------------------------------------------
function build_inertia()
%  Drag-sail pitch inertia matrices (kg m^2), about the CM, body axes.
%  (formerly drag_sail_inertias.m)
inertiaDB.angles_deg = (45:5:90).';
inertiaDB.I = NaN(3,3,numel(inertiaDB.angles_deg));

inertiaDB.I(:,:,inertiaDB.angles_deg==45) = diag([0.83 1.47 1.47]);
inertiaDB.I(:,:,inertiaDB.angles_deg==50) = diag([0.83 1.17 1.17]);
inertiaDB.I(:,:,inertiaDB.angles_deg==55) = diag([0.83 0.95 0.95]);
inertiaDB.I(:,:,inertiaDB.angles_deg==60) = diag([0.83 0.77 0.77]);
inertiaDB.I(:,:,inertiaDB.angles_deg==65) = diag([0.83 0.65 0.65]);
inertiaDB.I(:,:,inertiaDB.angles_deg==70) = diag([0.83 0.57 0.57]);
inertiaDB.I(:,:,inertiaDB.angles_deg==75) = diag([0.83 0.48 0.48]);
inertiaDB.I(:,:,inertiaDB.angles_deg==80) = diag([0.83 0.44 0.44]);
inertiaDB.I(:,:,inertiaDB.angles_deg==85) = diag([0.83 0.42 0.42]);
inertiaDB.I(:,:,inertiaDB.angles_deg==90) = diag([0.82 0.41 0.41]);

save('inertia_tensors.mat', 'inertiaDB');
fprintf('Wrote inertia_tensors.mat (%d apex angles)\n', numel(inertiaDB.angles_deg));
end


% -------------------------------------------------------------------------
function build_geometry(here)
%  Reference area / characteristic length per apex angle, from the .obj
%  geometries, using get_ref_geometry (ADBSat convention: A_ref = area/2).
phi_deg = 45:5:90;
A_ref = zeros(size(phi_deg));
L_ref = zeros(size(phi_deg));

for k = 1:numel(phi_deg)
    obj = fullfile(here, 'obj files', sprintf('%ddeg.obj', phi_deg(k)));
    if ~isfile(obj)
        warning('setup_databases:missingObj', 'Missing geometry: %s', obj);
        A_ref(k) = NaN; L_ref(k) = NaN;
        continue
    end
    [A_ref(k), L_ref(k)] = get_ref_geometry(obj);
end

ref_geom.phi_deg = phi_deg(:).';
ref_geom.A_ref   = A_ref(:).';
ref_geom.L_ref   = L_ref(:).';

save('ref_geometry.mat', 'ref_geom');
fprintf('Wrote ref_geometry.mat (%d apex angles)\n', numel(phi_deg));
end


% -------------------------------------------------------------------------
function rename_adbsat(here, resultsPath)
%  Copy timestamped ADBSat outputs into the expected folder/name layout.
%  (formerly postprocess.m)
phi_vec   = 45:5:90;
alpha_vec = 0.50:0.05:1.00;
alt_vec   = [350, 450, 650];
outBase   = fullfile(here, 'adbsat_processed');

for ai = 1:numel(alt_vec)
    alt    = alt_vec(ai);
    outDir = fullfile(outBase, sprintf('%dkm', alt));
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    for g = 1:numel(phi_vec)
        phi = phi_vec(g);
        for i = 1:numel(alpha_vec)
            a = alpha_vec(i);
            alpha_tag  = strrep(sprintf('alpha_%0.2f', a), '.', 'p');
            folderName = sprintf('%dkm_%ddeg_%s', alt, phi, alpha_tag);
            srcFolder  = fullfile(resultsPath, folderName);

            listing = dir(fullfile(srcFolder, '*.mat'));
            if isempty(listing)
                warning('setup_databases:noMat', 'No .mat in: %s', srcFolder);
                continue
            end
            srcFile   = fullfile(srcFolder, listing(1).name);
            accom_str = strrep(sprintf('%0.2f', a), '.', 'p');
            dstFile   = fullfile(outDir, sprintf('%ddeg_CLL_accom_%s.mat', phi, accom_str));
            copyfile(srcFile, dstFile);
            fprintf('Copied: %s\n', dstFile);
        end
    end
end
fprintf('\nDone. Files ready in: %s\n', outBase);
end
