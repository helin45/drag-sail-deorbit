function [A_ref, L_ref] = get_ref_geometry(obj_filepath)
% Reads an OBJ file and computes:
%   A_ref = total surface area / 2  (ADBSat convention)
%   L_ref = sqrt(A_ref)             (characteristic length)
%
% Usage:
%   [A_ref, L_ref] = get_ref_geometry(fullfile('obj files','45deg.obj'))

    fid = fopen(obj_filepath, 'r');
    if fid == -1
        error('Cannot open file: %s', obj_filepath);
    end

    verts = [];
    faces = [];

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if startsWith(line, 'v ')
            coords = sscanf(line(3:end), '%f %f %f')';
            verts  = [verts; coords];
        elseif startsWith(line, 'f ')
            % Handle both   f v1 v2 v3   and   f v1/vt1/vn1 ...
            tokens = strsplit(strtrim(line(3:end)));
            idx = zeros(1, length(tokens));
            for k = 1:length(tokens)
                parts  = strsplit(tokens{k}, '/');
                idx(k) = str2double(parts{1});
            end
            % Triangulate if face has more than 3 verts (fan method)
            for k = 2:length(idx)-1
                faces = [faces; idx(1), idx(k), idx(k+1)];
            end
        end
    end
    fclose(fid);

    % Compute total surface area
    total_area = 0;
    for f = 1:size(faces, 1)
        v1 = verts(faces(f,1), :);
        v2 = verts(faces(f,2), :);
        v3 = verts(faces(f,3), :);
        total_area = total_area + 0.5 * norm(cross(v2-v1, v3-v1));
    end

    A_ref = total_area / 2;
    L_ref = sqrt(A_ref);

    fprintf('File:   %s\n', obj_filepath);
    fprintf('A_ref = %.6f m^2\n', A_ref);
    fprintf('L_ref = %.6f m\n\n', L_ref);
end