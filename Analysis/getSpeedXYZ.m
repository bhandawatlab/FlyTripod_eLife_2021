function speed = getSpeedXYZ(xyz_mm, frame_rate)
    %getSpeedXYZ Determine speeds (mm/s) from a matrix of XYZ positions
    %(mm)
    xyz_diff = diff(xyz_mm, 1, 1);
    xdiff = xyz_diff(:, 1);
    ydiff = xyz_diff(:, 2);
    zdiff = xyz_diff(:, 3);
    dist_mm = sqrt(xdiff.^2 + ydiff.^2 + zdiff.^2);
    speed = dist_mm ./ (1/frame_rate);
end
