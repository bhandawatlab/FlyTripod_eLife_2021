function [data_xlim, data_ylim, data_zlim] =...
    getTrackingDataLimits(CTr_XYZ, FTi_XYZ, TiTa_XYZ, Ta_XYZ)
    %GETTRACKINGDATALIMITS Get the data's XYZ limits.
    data_xlim = [NaN NaN];
    data_ylim = [NaN NaN];
    data_zlim = [NaN NaN];
    
    xyz_data_sets = {CTr_XYZ, FTi_XYZ, TiTa_XYZ, Ta_XYZ};
    for n = 1:length(xyz_data_sets)
        xyz_data = xyz_data_sets{n};
        
        % Update the limits
        data_xlim(1) = min(data_xlim(1), min(xyz_data(:, 1)));
        data_xlim(2) = max(data_xlim(2), max(xyz_data(:, 1)));
        data_ylim(1) = min(data_ylim(1), min(xyz_data(:, 2)));
        data_ylim(2) = max(data_ylim(2), max(xyz_data(:, 2)));
        data_zlim(1) = min(data_zlim(1), min(xyz_data(:, 3)));
        data_zlim(2) = max(data_zlim(2), max(xyz_data(:, 3)));
    end
end

