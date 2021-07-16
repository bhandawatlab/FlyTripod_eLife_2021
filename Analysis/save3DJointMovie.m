function save3DJointMovie(tracking_data_file, tracking_data)
    %save3DJointMovie Save the 3D joint positions as a movie.
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
    
    % Save a movie of the joint positions
    % (Reference: https://www.mathworks.com/help/matlab/ref/getframe.html)
    % https://www.mathworks.com/help/matlab/ref/movie.html
    
    h = figure;
%     axis tight manual
    axis tight
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    
    % Preallocate an array M to store the movie frames.
    frame_count = size(tracking_data.Anterior_XYZ, 1);
%     frame_count = 10;
    M(frame_count) = struct('cdata',[],'colormap',[]);
    
    % Capture each plot of function X as an individual frame and store them in M.
    % Set the 'Visible' property of the figure object to 'off'.
    h.Visible = 'off';
    bad_frames = tracking_data.bad_frames;
    for frame_index = 1:frame_count
        if ~bad_frames(frame_index)
            % TODO: Center all points around the first joint (coxa) so that
            % it stays fixed in the movie.
            
            %% Plot the leg
            leg_x = [];
            leg_y = [];
            leg_z = [];
            for joint_index = 1:joint_count
                joint_name = joint_names{joint_index};
                joint_xyz_array = tracking_data.(joint_name);
                joint_xyz_point = joint_xyz_array(frame_index,:);
                leg_x = [leg_x; joint_xyz_point(1)];
                leg_y = [leg_y; joint_xyz_point(2)];
                leg_z = [leg_z; joint_xyz_point(3)];

                % Plot the joint
                joint_color = joint_colors{joint_index};
                scatter3(joint_xyz_point(1), joint_xyz_point(2), joint_xyz_point(3), 'MarkerFaceColor', joint_color, 'MarkerEdgeColor', 'k');
                hold on
            end

            % Plot the leg
            plot3(leg_x, leg_y, leg_z, 'k');

%             %% Plot the anterior-posterior line as a dashed line
%             anterior_xyz_array = tracking_data.('Anterior_XYZ');
%             anterior_xyz_point = anterior_xyz_array(frame_index, :);
%             posterior_xyz_array = tracking_data.('Posterior_XYZ');
%             posterior_xyz_point = posterior_xyz_array(frame_index, :);
%             ap_x = [anterior_xyz_point(1); posterior_xyz_point(1)];
%             ap_y = [anterior_xyz_point(2); posterior_xyz_point(2)];
%             ap_z = [anterior_xyz_point(3); posterior_xyz_point(3)];
%             plot3(ap_x, ap_y, ap_z, 'k--');
% 
%             %% Plot the anterior point as a filled white circle
%             scatter3(anterior_xyz_point(1), anterior_xyz_point(2), anterior_xyz_point(3), 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
%             hold off

            %% Set the aspect ratio as equal
            daspect([1 1 1])

            %% Update and grab the figure as a frame
            drawnow
            M(frame_index) = getframe;
        else
            M(frame_index) = getframe;
        end
    end
    
    % Set the 'Visible' property of the figure to 'on'. Play the movie once at 6 frames per second.
    h.Visible = 'on';
    
    movie(M,1,6);
    %%%%%%%
%     figure;
%     
%     %% Plot the anterior-posterior points
%     anterior_xyz = tracking_data.('Anterior_XYZ');
%     posterior_xyz = tracking_data.('Posterior_XYZ');
%     plot3(anterior_xyz(:,1),anterior_xyz(:,2),anterior_xyz(:,3), 'k.')
%     hold on
%     plot3(posterior_xyz(:,1),posterior_xyz(:,2),posterior_xyz(:,3), 'ko')
% 
%     %% Plot the joints
%     for joint_index = 1:length(joint_names)
%         joint_name = joint_names{joint_index};
%         joint_color = joint_colors{joint_index};
%         joint_xyz = tracking_data.(joint_name);
%         
%         % Plot
%         plot3(joint_xyz(:,1),joint_xyz(:,2),joint_xyz(:,3), 'Color', joint_color)
%     end
%     hold off
%     daspect([1 1 1])
end

