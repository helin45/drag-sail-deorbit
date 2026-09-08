%% =========================================================================
%  POST-PROCESSING: Rename ADBSat outputs -> attitude dynamics model format
%  Renames:  <timestamp>.mat  ->  <angle>deg_CLL_accom_<accom>.mat
%  Organised into subfolders per altitude in your working directory
% =========================================================================
clear; clc;

resultsPath = '/Users/helintaha/Library/CloudStorage/OneDrive-TheUniversityofManchester/Dissertation/ADBSat-master/inou/results';

phi_vec   = 45:5:90;
alpha_vec = 0.50:0.05:1.00;
alt_vec   = [350, 450, 650];

outBase = fullfile(pwd, 'adbsat_processed');

for ai = 1:length(alt_vec)
    alt    = alt_vec(ai);
    outDir = fullfile(outBase, sprintf('%dkm', alt));
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    for g = 1:length(phi_vec)
        phi     = phi_vec(g);
        modName = sprintf('%ddeg', phi);

        for i = 1:length(alpha_vec)
            a = alpha_vec(i);

            alpha_tag  = strrep(sprintf('alpha_%0.2f', a), '.', 'p');
            folderName = sprintf('%dkm_%s_%s', alt, modName, alpha_tag);
            srcFolder  = fullfile(resultsPath, folderName);

            % Find the .mat file (timestamped name)
            listing = dir(fullfile(srcFolder, '*.mat'));
            if isempty(listing)
                warning('No .mat found in: %s', srcFolder);
                continue
            end

            srcFile = fullfile(srcFolder, listing(1).name);

            % Destination filename matches what attitude model expects
            accom_str = strrep(sprintf('%0.2f', a), '.', 'p');
            dstFile   = fullfile(outDir, ...
                        sprintf('%ddeg_CLL_accom_%s.mat', phi, accom_str));

            copyfile(srcFile, dstFile);
            fprintf('Copied: %s\n', dstFile);
        end
    end
end

fprintf('\nDone. Files ready in: %s\n', outBase);