function root = project_root()
%PROJECT_ROOT  Absolute path to the repository root folder.
%
%   root = PROJECT_ROOT() returns the folder that contains run_sweep.m,
%   run_analyses.m, make_figures.m and the adbsat_processed/ data tree.
%
%   This file lives in <root>/lib, so the root is one level up. Every script
%   under lib/analyses and lib/figures uses this instead of a hard-coded
%   absolute path, so the repo works wherever it is cloned.

root = fileparts(fileparts(mfilename('fullpath')));
end
