function plot3DJointPositions(tracking_data)
    %PLOT3DJOINTPOSITIONS Plot the 3D joint positions.
    % Colors for the five joints
    joint_colors = {'#A2142F' '#77AC30' '#7E2F8E' '#0072BD' '#D95319'};
    
    % Set up a list (and order) of joints to plot
    joint_names = {...
        'CTr_R_Pro_XYZ'
        'FTi_R_Pro_XYZ'
        'TiTa_R_Pro_XYZ'
        'Ta_XYZ'
        };
    
    % Save a movie of the joint positions
    % (Reference: https://www.mathworks.com/help/matlab/ref/getframe.html)
    % https://www.mathworks.com/help/matlab/ref/movie.html
    
    figure;
    
    %% Plot the anterior-posterior line
    anterior_xyz = tracking_data.('Anterior_XYZ');
    posterior_xyz = tracking_data.('Posterior_XYZ');
    ap_line = [anterior_xyz; posterior_xyz
    plot3(CTr_R_Pro_XYZ_mm(:,1),CTr_R_Pro_XYZ_mm(:,2),CTr_R_Pro_XYZ_mm(:,3))
end

