function tracking_data = updateFly3DTrackingData(tracking_data_file, mm_per_pixel, frame_rate, max_speed, pcutoff)
    %updateFly3DTrackingData Update the fly's tracking data file with the
    %3D-reconstructed points
    % Set up a list of pairs of points to reconstruct (Bottom view in the
    % first column, top view in the 2nd column)
    point_pairs = {...
        'Anterior_Bottom' 'Anterior_Top'
        'CTr_R_Pro_Bottom' 'CTr_R_Pro_Top'
        'FTi_R_Pro_Bottom' 'FTi_R_Pro_Top'
        'Posterior_Bottom' 'Posterior_Top'
        'Ta_Bottom' 'Ta_Top'
        'TiTa_R_Pro_Bottom' 'TiTa_R_Pro_Top'
        };
    
    % Load the tracking data
    tracking_data = load(tracking_data_file, '-mat');
    
    % Keep track of all bad frames
    frame_count = length(tracking_data.('Anterior_Bottom').x);
    bad_frame_indices = zeros(frame_count-1, 1);
    
    % Loop through all pairs of points, 3D-reconstruct, and store the result
    for pair_index = 1:size(point_pairs, 1)
        % Field name for bottom view
        bottom_view_name = point_pairs{pair_index, 1};
        
        % Field name for top view
        top_view_name = point_pairs{pair_index, 2};
        
        %% Access i,j pixel positions for both sets of points
        bottom_ij_data = tracking_data.(bottom_view_name);
        bottom_ij = [bottom_ij_data.x' bottom_ij_data.y'];
        top_ij_data = tracking_data.(top_view_name);
        top_ij = [top_ij_data.x' top_ij_data.y'];
        
        % Set to NaN all tracks with low likelihood
        bottom_ij_lh = bottom_ij_data.likelihood;
        bottom_ij(bottom_ij_lh <= pcutoff, 1:2) = NaN;
        top_ij_lh = top_ij_data.likelihood;
        top_ij(top_ij_lh <= pcutoff, 1:2) = NaN;
        
        %% Reconstruct into world coordinates (mm)
        xyz_mm = Mirror3DReconstruction(top_ij, bottom_ij, mm_per_pixel);

        %% Ignore all values with speeds > 30mm/s (tracking errors)
        % Get the speed (mm/s)
        xyz_speed = getSpeedXYZ(xyz_mm, frame_rate);
        
        % Remove the first index to compare directly with speed
        xyz_mm = xyz_mm(2:end,:);
        
        % Set high-speed values to NaN
        xyz_mm(xyz_speed >= max_speed) = NaN;
        
        %% Add the XYZ points to the tracking data
        % Determine the field name for storing 3D-reconstructed points
        split_top_view_name = split(top_view_name, 'Top');
        reconstructed_name = [split_top_view_name{1} 'XYZ'];
        
        % Store the XYZ field
        tracking_data.(reconstructed_name) = xyz_mm;
        
        %% Update the bad frame indices
        nan_positions = isnan(xyz_mm);
        nan_logical = logical(sum(nan_positions, 2));
        bad_frame_indices = or(bad_frame_indices, nan_logical);
    end
    
    %% Add the bad frames to the tracking data
    tracking_data.bad_frames = bad_frame_indices;
    
    %% Save the updated tracking file
    save(tracking_data_file, '-struct', 'tracking_data');
end
