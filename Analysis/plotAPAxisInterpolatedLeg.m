function plotAPAxisInterpolatedLeg(ap_axis_points, tracking_data_file)
    %plotAPAxisInterpolatedLeg Plot the leg relative to the
    %anterior-posterior body axis. Save an animated *.GIF of this plot.
    % Get the leg points
    ThC_U = ap_axis_points(1).ThC_U;
    CTr_U = ap_axis_points(1).CTr_U;
    FTi_U = ap_axis_points(1).FTi_U;
    TiTa_U = ap_axis_points(1).TiTa_U;
    Ta_U = ap_axis_points(1).Ta_U;
    
    % Determine the output *.GIF filename
    split_filename = split(tracking_data_file, '.mat');
    gif_filename = [split_filename{1} '_AP-axis-Interp-Test.gif'];
    gif_xz_filename = [split_filename{1} '_AP-axis-XZ-Interp-Test.gif'];
    
    % Colors  and names for the four labels
    label_colors = {'#A2142F' '#77AC30' '#7E2F8E' '#0072BD'};
    label_names = {...
        'ThC_R_Pro_XYZ'
        'CTr_R_Pro_XYZ'
        'FTi_R_Pro_XYZ'
        'TiTa_R_Pro_XYZ'
        'Ta_XYZ'
        };
        
    %% Determine the plot limits
    [data_xlim, data_ylim, data_zlim] =...
        getTrackingDataLimits(CTr_U, FTi_U, TiTa_U, Ta_U);
    
    %% Plot
    % Initialize the figure
    main_figure = figure;
    axis tight manual
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    set(gcf,'renderer','painters')
    point_count = size(CTr_U, 1);
    for n=1:point_count
        % Plot the leg
        leg_x = [ThC_U(1,1); CTr_U(n,1); FTi_U(n,1); TiTa_U(n,1); Ta_U(n,1)];
        leg_y = [ThC_U(1,2); CTr_U(n,2); FTi_U(n,2); TiTa_U(n,2); Ta_U(n,2)];
        leg_z = [ThC_U(1,3); CTr_U(n,3); FTi_U(n,3); TiTa_U(n,3); Ta_U(n,3)];
        h = plot3(leg_x, leg_y, leg_z, 'k');
        hold on

        % Turn off legends for the leg
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';

        % Plot the joints and tarsus tip
        scatter3(ThC_U(1,1), ThC_U(1,2), ThC_U(1,3), 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'k');
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
end
