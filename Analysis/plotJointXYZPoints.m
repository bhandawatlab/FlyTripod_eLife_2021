function plotJointXYZPoints(tracking_data)
    %PLOTJOINTXYZPOINTS Plot joint positions.

    % Set up a list (and order) of points to plot
    field_names = {...
        'Anterior_XYZ'
        'Posterior_XYZ'
        'CTr_R_Pro_XYZ'
        'FTi_R_Pro_XYZ'
        'TiTa_R_Pro_XYZ'
        'Ta_XYZ'
        };
    field_count = length(field_names);

    % Access the good frame indices
    good_frame_indices = ~tracking_data.bad_frames;
    
    % Generate the XY (front view) and YZ (side view) plots
    for field_index = 1:field_count
        % Access the points for this field
        field_name = field_names{field_index};
        xyz_array = tracking_data.(field_name);
        xyz_array = xyz_array(good_frame_indices, :); % Use good frames only

        % Plot the field's points
        f = figure;
        plot_layout = tiledlayout(f, 1, 2);
        
        % XY
        xy_axis = nexttile;
        scatter(xy_axis, xyz_array(:, 1), xyz_array(:, 2), 'filled');
        title(xy_axis, 'XY (Front view)')
        
        % YZ
        yz_axis = nexttile;
        scatter(yz_axis, xyz_array(:, 2), xyz_array(:, 3), 'filled');
        title(yz_axis, 'YZ (Side view)')
        title(plot_layout, field_name, 'Interpreter', 'none')
    end
    
%     %% Plot the anterior-posterior points
%     anterior_xyz = tracking_data.('Anterior_XYZ');
%     posterior_xyz = tracking_data.('Posterior_XYZ');
%     plot3(anterior_xyz(:,1),anterior_xyz(:,2),anterior_xyz(:,3), 'k.')
%     plot3(posterior_xyz(:,1),posterior_xyz(:,2),posterior_xyz(:,3), 'ko')
%     hold off
end
