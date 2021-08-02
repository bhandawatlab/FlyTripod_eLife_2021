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
    bad_frame_logical = zeros(frame_count-1, 1);
    
    % Loop through all pairs of points, 3D-reconstruct, and store the result
    point_set_count = size(point_pairs, 1);
    xyz_label_names = cell(point_set_count, 1);
    for pair_index = 1:point_set_count
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
        xyz_label_names{pair_index} = reconstructed_name;
        
        % Store the XYZ field
        tracking_data.(reconstructed_name) = xyz_mm;
        
        %% Update the bad frame indices
        nan_positions = isnan(xyz_mm);
        nan_logical = logical(sum(nan_positions, 2));
        bad_frame_logical = or(bad_frame_logical, nan_logical);
    end
    
    %% Add the bad frames to the tracking data
    tracking_data.bad_frames = bad_frame_logical;
    
    %% Interpolate over the bad frames
    % Get the good frame start and end
    good_frame_logical = ~bad_frame_logical;
    good_frame_indices = find(good_frame_logical);
    start_good_frame = good_frame_indices(1);
    end_good_frame = good_frame_indices(end);
    
    % Ignore NaNs at the start and end
    good_frame_logical = good_frame_logical(start_good_frame:end_good_frame);
    bad_frame_logical = bad_frame_logical(start_good_frame:end_good_frame);
    
    % Loop through all points
    for xyz_label_index = 1:point_set_count
        % Get the XYZ points
        xyz_label_name = xyz_label_names{xyz_label_index};
        xyz_mm = tracking_data.(xyz_label_name);
        
        % Ignore NaNs at the start and end
        xyz_mm = xyz_mm(start_good_frame:end_good_frame, :);
        
        % Interpolate over NaNs (1D)
        indices = 1:size(xyz_mm, 1);
        sample_indices = indices(good_frame_logical);
        sample_values = xyz_mm(good_frame_logical, :);
        query_indices = indices(bad_frame_logical);
        x_interp = interp1(sample_indices, sample_values(:, 1), query_indices);
        y_interp = interp1(sample_indices, sample_values(:, 2), query_indices);
        z_interp = interp1(sample_indices, sample_values(:, 3), query_indices);
        xyz_mm(query_indices, :) = [x_interp' y_interp' z_interp'];
        
        % Add the interpolated points to the output structure
        interpolated_name = [xyz_label_name '_Interp'];
        tracking_data.(interpolated_name) = xyz_mm;
    end
    
    %% Save the updated tracking file
    save(tracking_data_file, '-struct', 'tracking_data');
end
