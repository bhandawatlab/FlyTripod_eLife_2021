function plotLegAtAnteriorPosteriorAxis(tracking_data_file, tracking_data, coxa_length)
    %plotLegAtAnteriorPosteriorAxis Plot the leg relative to the
    %anterior-posterior body axis. Save an animated *.GIF of this plot.
    
    % Determine the output *.GIF filename
    split_filename = split(tracking_data_file, '.mat');
    gif_filename = [split_filename{1} '_AP-axis.gif'];
    gif_xz_filename = [split_filename{1} '_AP-axis-XZ.gif'];
    
    % Colors  and names for the four labels
    label_colors = {'#A2142F' '#77AC30' '#7E2F8E' '#0072BD'};
    label_names = {...
        'CTr_R_Pro_XYZ'
        'FTi_R_Pro_XYZ'
        'TiTa_R_Pro_XYZ'
        'Ta_XYZ'
        };
    
    % Access all good points
    good_frame_indices = find(~tracking_data.bad_frames);
    Anterior_XYZ = tracking_data.Anterior_XYZ(good_frame_indices, :);
    Posterior_XYZ = tracking_data.Posterior_XYZ(good_frame_indices, :);
    CTr_XYZ = tracking_data.CTr_R_Pro_XYZ(good_frame_indices, :);
    FTi_XYZ = tracking_data.FTi_R_Pro_XYZ(good_frame_indices, :);
    TiTa_XYZ = tracking_data.TiTa_R_Pro_XYZ(good_frame_indices, :);
    Ta_XYZ = tracking_data.Ta_XYZ(good_frame_indices, :);
    
    % Define the posterior-anterior vector relative to the posterior point
    PA_vec = Anterior_XYZ - Posterior_XYZ;
    
    % Define the position vectors relative to the posterior
    % point. The origin is now at the posterior point.
    CTr_XYZ_PA = CTr_XYZ - Posterior_XYZ;
    FTi_XYZ_PA = FTi_XYZ - Posterior_XYZ;
    TiTa_XYZ_PA = TiTa_XYZ - Posterior_XYZ;
    Ta_XYZ_PA = Ta_XYZ - Posterior_XYZ;
    
    % Define new normalized bases for each point relative to the AP-axis
    % using the double cross product method
    % u1 is the AP-axis
    u1 = PA_vec ./ vecnorm(PA_vec, 2, 2);
    point_count = size(CTr_XYZ, 1);
    y_axis = repmat([0 1 0], point_count, 1);
    
    % u2 is the vector orthogonal to the standard Y (DV) and u1
    % axes
    u2 = cross(u1, y_axis, 2);
    u2 = u2 / norm(u2);
    
    % u3 is the vector orthogonal to the u1 and u2 axes
    u3 = cross(u1, u2, 2);
    u3 = u3 / norm(u3);
    
    % Plot the axes
%     plotAxes(u1(1, :), u2(1, :), u3(1, :));

    % Loop for transforming the CTr positions to the AP vector bases
    CTr_U = NaN(point_count, 3);
    FTi_U = NaN(point_count, 3);
    TiTa_U = NaN(point_count, 3);
    Ta_U = NaN(point_count, 3);
    for n=1:point_count
        % Get the basis vectors
        u1_vec = u1(n, :);
        u2_vec = u2(n, :);
        u3_vec = u3(n, :);
                
        % Set up the transition matrix between bases
        U = [u1_vec' u2_vec' u3_vec'];
        U_inv = inv(U);

        % Transform the CTr positions to the AP vector basis
        CTr_XYZ_PA_vec = CTr_XYZ_PA(n, :);
        CTr_U_vec = U_inv * CTr_XYZ_PA_vec';
        CTr_U(n, :) = CTr_U_vec;
        
        % Transform the FTi positions to the AP vector basis
        FTi_XYZ_PA_vec = FTi_XYZ_PA(n, :);
        FTi_U_vec = U_inv * FTi_XYZ_PA_vec';
        FTi_U(n, :) = FTi_U_vec;
        
        % Transform the TiTa positions to the AP vector basis
        TiTa_XYZ_PA_vec = TiTa_XYZ_PA(n, :);
        TiTa_U_vec = U_inv * TiTa_XYZ_PA_vec';
        TiTa_U(n, :) = TiTa_U_vec;
        
        % Transform the Ta positions to the AP vector basis
        Ta_XYZ_PA_vec = Ta_XYZ_PA(n, :);
        Ta_U_vec = U_inv * Ta_XYZ_PA_vec';
        Ta_U(n, :) = Ta_U_vec;
    end
    
    % Determine the plot limits
    [data_xlim, data_ylim, data_zlim] =...
        getTrackingDataLimits(CTr_U, FTi_U, TiTa_U, Ta_U);
    
    %% Plot
    % Initialize the figure
    main_figure = figure;
    axis tight manual
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    set(gcf,'renderer','painters')
    for n=1:point_count        
        % Plot the leg
        leg_x = [CTr_U(n,1); FTi_U(n,1); TiTa_U(n,1); Ta_U(n,1)];
        leg_y = [CTr_U(n,2); FTi_U(n,2); TiTa_U(n,2); Ta_U(n,2)];
        leg_z = [CTr_U(n,3); FTi_U(n,3); TiTa_U(n,3); Ta_U(n,3)];
        h = plot3(leg_x, leg_y, leg_z, 'k');
        hold on
        
        % Turn off legends for the leg
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        % Plot the joints and tarsus tip
        scatter3(CTr_U(n,1), CTr_U(n,2), CTr_U(n,3), 'MarkerFaceColor', label_colors{1}, 'MarkerEdgeColor', 'k');
        scatter3(FTi_U(n,1), FTi_U(n,2), FTi_U(n,3), 'MarkerFaceColor', label_colors{2}, 'MarkerEdgeColor', 'k');
        scatter3(TiTa_U(n,1), TiTa_U(n,2), TiTa_U(n,3), 'MarkerFaceColor', label_colors{3}, 'MarkerEdgeColor', 'k');
        scatter3(Ta_U(n,1), Ta_U(n,2), Ta_U(n,3), 'MarkerFaceColor', label_colors{4}, 'MarkerEdgeColor', 'k');
        legend(label_names, 'Location', 'southeastoutside', 'Interpreter', 'none')
        hold off
        
        % Set the plot limits
        xlim(data_xlim)
        ylim(data_ylim)
        zlim(data_zlim)
        daspect([1 1 1])  % Set the aspect ratio as equal
        xlabel('X (AP)')
        ylabel('Y (LR)')
        zlabel('Z (DV)')
        
        %% Save the 3D view as a *.GIF (https://www.mathworks.com/help/matlab/ref/imwrite.html#btv452g-1)
        % Capture the plot as an image 
        drawnow
        current_frame = getframe(main_figure); 
        im = frame2im(current_frame); 
        [imind,cm] = rgb2ind(im,256);

        if n == 1
            imwrite(imind, cm, gif_filename,'gif','LoopCount',Inf,'DelayTime',1/30);
        else
            imwrite(imind, cm, gif_filename,'gif','WriteMode','append','DelayTime',1/30);
        end
        
        %% Save the 2D XZ view as a *.GIF (https://www.mathworks.com/help/matlab/ref/imwrite.html#btv452g-1)
        % Rotate the view to align it into the XZ axis (Z pointing
        % vertically)
        % (https://www.mathworks.com/help/matlab/ref/view.html#d123e1472755)
        view(0,0);
        
        % Capture the plot as an image 
        drawnow
        current_frame = getframe(main_figure); 
        im = frame2im(current_frame); 
        [imind,cm] = rgb2ind(im,256);

        if n == 1
            imwrite(imind, cm, gif_xz_filename,'gif','LoopCount',Inf,'DelayTime',1/30);
        else
            imwrite(imind, cm, gif_xz_filename,'gif','WriteMode','append','DelayTime',1/30);
        end
    end
    
%     % Plot the transformed CTr vectors to check
%     figure;
%     scatter3(CTr_U(:, 1), CTr_U(:, 2), CTr_U(:, 3), 'filled');
%     hold on
%     scatter3(0,0,0, 'filled', 'r');
%     hold off
end
