% With the fly facing the right side of the top view:
% X = A/P axis
% Y = D/V axis
% Z = L/R axis (depth) = Must be calculated from the mirror plane
% [Top view = XY plane, Mirror view = XZ plane]
% Calculate the 3D reconstruction given a set of points XY in the top
% view, and XZ in the bottom 45-degree mirror view.


function xyz_mm = Mirror3DReconstruction(ij_top, ij_bottom, mm_per_pixel)
    % Inputs are the top and bottom view coordinates in pixels, as well as
    % the pixel spacing in mm.
    
    % The top view is our ground truth. The XY position can be determined
    % from the pixel spacing:
    xy_mm = ij_top .* mm_per_pixel;
    
    % The depth is determined from the mirror reference frame and the
    % mirrored Y:
    y2 = ij_bottom(:, 2);
    r_t = 170;  % Distance (mm) from the lens to the back of the chamber
    mirror_theta = pi/4;  % Mirror angle
    z_num = y2 - (r_t/mm_per_pixel)*tan((2*mirror_theta) - (pi/2));
    z_den = tan(mirror_theta) - tan(2*mirror_theta - (pi/2));
    z_mm = (z_num ./ z_den) .* mm_per_pixel;
    
    % Format the output array
    xyz_mm = [xy_mm z_mm];
end
