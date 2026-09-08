% drag_sail_inertias.m
% Drag sail inertia matrices (kg*m^2), about CM, aligned with body axes

inertiaDB.angles_deg = (45:5:90).';                 % 10x1
inertiaDB.I = NaN(3,3,numel(inertiaDB.angles_deg)); % 3x3x10

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

% Save as a .mat that SAM_from_ADBSat loads
save("inertia_tensors.mat","inertiaDB");