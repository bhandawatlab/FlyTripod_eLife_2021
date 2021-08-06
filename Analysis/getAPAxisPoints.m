function ap_axis_points = getAPAxisPoints(tracking_data)
    %getAPAxisPoints Get the leg points relative to the anterior-posterior
    %body axis.
    % Set up the output structure of transformed leg points
    ap_axis_points = struct();
    
    % Access all interpolated points
    Anterior_XYZ = tracking_data.Anterior_XYZ_Interp;
    Posterior_XYZ = tracking_data.Posterior_XYZ_Interp;
    CTr_XYZ = tracking_data.CTr_R_Pro_XYZ_Interp;
    FTi_XYZ = tracking_data.FTi_R_Pro_XYZ_Interp;
    TiTa_XYZ = tracking_data.TiTa_R_Pro_XYZ_Interp;
    Ta_XYZ = tracking_data.Ta_XYZ_Interp;

    % Define the posterior-anterior vector relative to the posterior point
    PA_XYZ = Anterior_XYZ - Posterior_XYZ;
    
    % Define the position vectors relative to the posterior
    % point. The origin is now at the posterior point.
    CTr_XYZ_PA = CTr_XYZ - Posterior_XYZ;
    FTi_XYZ_PA = FTi_XYZ - Posterior_XYZ;
    TiTa_XYZ_PA = TiTa_XYZ - Posterior_XYZ;
    Ta_XYZ_PA = Ta_XYZ - Posterior_XYZ;

    %% Loop for transforming the CTr positions to the AP vector bases
    point_count = size(PA_XYZ, 1);
    CTr_U = NaN(point_count, 3);
    FTi_U = NaN(point_count, 3);
    TiTa_U = NaN(point_count, 3);
    Ta_U = NaN(point_count, 3);
    for n=1:point_count
        % Get the anterior-posterior axis
        PA_XYZ_vec = PA_XYZ(n, :);

        % Define new normalized bases for each point relative to the
        % AP-axis using the double cross product method
        % u1 is the AP-axis
        u1 = PA_XYZ_vec / norm(PA_XYZ_vec);
        y_axis = [0 1 0];

        % u2 is the vector orthogonal to the standard Y (DV) and u1
        % axes
        u2 = cross(u1, y_axis);
        u2 = u2 / norm(u2);

        % u3 is the vector orthogonal to the u1 and u2 axes
        u3 = cross(u1, u2);
        u3 = u3 / norm(u3);

        % Plot the axes
        %plotAxes(u1(1, :), u2(1, :), u3(1, :));

        % Set up the transition matrix between bases
        U = [u1' u2' u3'];
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
    
    %% Find the ThC joint by fitting a sphere to the CTr point cloud
    CTr_U_point_cloud = pointCloud(CTr_U);
    max_distance = 500;
    sphere_model = pcfitsphere(CTr_U_point_cloud, max_distance);
    ThC_U = sphere_model.Center;
    coxa_length_um = sphere_model.Radius * 1000;
    sprintf('Coxa length: %.3f µm\n', coxa_length_um);
    
    %     % Plot the sphere
    %     figure;
    %     plot(sphere_model);
    %     alpha 0
    %     hold on
    %     scatter3(CTr_U(:,1), CTr_U(:,2), CTr_U(:,3), 'filled')
    %     hold off
    
    %% Store the transformed leg points
    ap_axis_points(1).ThC_U = ThC_U;
    ap_axis_points(1).CTr_U = CTr_U;
    ap_axis_points(1).FTi_U = FTi_U;
    ap_axis_points(1).TiTa_U = TiTa_U;
    ap_axis_points(1).Ta_U = Ta_U;
end
