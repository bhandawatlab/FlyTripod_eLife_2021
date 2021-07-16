function plot3DJointPositions(tracking_data)
    %PLOT3DJOINTPOSITIONS Plot the 3D joint positions.
%     % Colors for the five labels
%     joint_colors = {'#A2142F' '#77AC30' '#7E2F8E' '#0072BD' '#D95319'};
    % Colors for the four labels
    joint_colors = {'#A2142F' '#77AC30' '#7E2F8E' '#0072BD'};

    % Set up a list (and order) of joints to plot
    joint_names = {...
        'CTr_R_Pro_XYZ'
        'FTi_R_Pro_XYZ'
        'TiTa_R_Pro_XYZ'
        'Ta_XYZ'
        };
    
    figure;
    
    %% Plot the anterior-posterior points
    anterior_xyz = tracking_data.('Anterior_XYZ');
    posterior_xyz = tracking_data.('Posterior_XYZ');
    plot3(anterior_xyz(:,1),anterior_xyz(:,2),anterior_xyz(:,3), 'k.')
    hold on
    plot3(posterior_xyz(:,1),posterior_xyz(:,2),posterior_xyz(:,3), 'ko')

    %% Plot the joints
    for joint_index = 1:length(joint_names)
        joint_name = joint_names{joint_index};
        joint_color = joint_colors{joint_index};
        joint_xyz = tracking_data.(joint_name);
        
        % Plot
        plot3(joint_xyz(:,1),joint_xyz(:,2),joint_xyz(:,3), 'Color', joint_color)
    end
    hold off
    daspect([1 1 1])
end

