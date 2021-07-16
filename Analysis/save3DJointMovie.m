function save3DJointMovie(tracking_data_file, tracking_data)
    %save3DJointMovie Save the 3D joint positions as a *.GIF.
    % Colors for the four labels
    joint_colors = {'#A2142F' '#77AC30' '#7E2F8E' '#0072BD'};

    % Set up a list (and order) of joints to plot
    joint_names = {...
        'CTr_R_Pro_XYZ'
        'FTi_R_Pro_XYZ'
        'TiTa_R_Pro_XYZ'
        'Ta_XYZ'
        };
    joint_count = length(joint_names);
    
    % Keep track of the XYZ limits across joints
    xlim_joints = [NaN NaN];
    ylim_joints = [NaN NaN];
    zlim_joints = [NaN NaN];
    
    % Initialize the figure
    main_figure = figure;
    axis tight manual
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    set(gcf,'renderer','painters')
    
    % Preallocate an array M to store the movie frames.
    good_frame_indices = find(~tracking_data.bad_frames);
    good_frame_count = length(good_frame_indices);
    
    % Capture each plot of function X as an individual frame and store them in M.
    % Set the 'Visible' property of the figure object to 'off'.
%     h.Visible = 'off';
    
    %% First determine the plot limits
    for n = 1:good_frame_count
        frame_index = good_frame_indices(n);
        %% Plot the leg
        for joint_index = 1:joint_count
            % Access the joint position
            joint_name = joint_names{joint_index};
            joint_xyz_array = tracking_data.(joint_name);
            joint_xyz_point = joint_xyz_array(frame_index,:);

            % Update the limits
            xlim_joints(1) = min(xlim_joints(1), joint_xyz_point(1));
            xlim_joints(2) = max(xlim_joints(2), joint_xyz_point(1));
            ylim_joints(1) = min(ylim_joints(1), joint_xyz_point(2));
            ylim_joints(2) = max(ylim_joints(2), joint_xyz_point(2));
            zlim_joints(1) = min(zlim_joints(1), joint_xyz_point(3));
            zlim_joints(2) = max(zlim_joints(2), joint_xyz_point(3));
        end
    end
    
    %% Save a *.GIF of the leg positions
    % Determine the filename
    split_filename = split(tracking_data_file, '.mat');
    gif_filename = [split_filename{1} '.gif'];
    
    % Generate the plots
    for n = 1:good_frame_count
        frame_index = good_frame_indices(n);
        
        %% Plot the leg
        leg_x = [];
        leg_y = [];
        leg_z = [];
        for joint_index = 1:joint_count
            % Access the joint position
            joint_name = joint_names{joint_index};
            joint_xyz_array = tracking_data.(joint_name);
            joint_xyz_point = joint_xyz_array(frame_index,:);

            % Update the leg points
            leg_x = [leg_x; joint_xyz_point(1)];
            leg_y = [leg_y; joint_xyz_point(2)];
            leg_z = [leg_z; joint_xyz_point(3)];

            % Plot the joint
            joint_color = joint_colors{joint_index};
            scatter3(joint_xyz_point(1), joint_xyz_point(2), joint_xyz_point(3), 'MarkerFaceColor', joint_color, 'MarkerEdgeColor', 'k');
            xlim(xlim_joints)
            ylim(ylim_joints)
            zlim(zlim_joints)
            hold on
        end
        legend(joint_names, 'Location', 'southeast', 'Interpreter', 'none')

        % Plot the leg
        h = plot3(leg_x, leg_y, leg_z, 'k');
        
        % Turn off the legend
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        % Set the plot limits
        xlim(xlim_joints)
        ylim(ylim_joints)
        zlim(zlim_joints)
        hold off

        % Set the aspect ratio as equal
        daspect([1 1 1])
        
        % Flip the Z-axis
        set(gca, 'ZDir','reverse')
        set(gca, 'YDir','reverse')

        % Capture the plot as an image 
        drawnow
        current_frame = getframe(main_figure); 
        im = frame2im(current_frame); 
        [imind,cm] = rgb2ind(im,256);

        % Save as a *.GIF (https://www.mathworks.com/help/matlab/ref/imwrite.html#btv452g-1)
        if n == 1
            imwrite(imind, cm, gif_filename,'gif','LoopCount',Inf,'DelayTime',1/30);
        else
            imwrite(imind, cm, gif_filename,'gif','WriteMode','append','DelayTime',1/30);
        end
    end
    
%     %% Plot the anterior-posterior points
%     anterior_xyz = tracking_data.('Anterior_XYZ');
%     posterior_xyz = tracking_data.('Posterior_XYZ');
%     plot3(anterior_xyz(:,1),anterior_xyz(:,2),anterior_xyz(:,3), 'k.')
%     plot3(posterior_xyz(:,1),posterior_xyz(:,2),posterior_xyz(:,3), 'ko')
%     hold off
end
