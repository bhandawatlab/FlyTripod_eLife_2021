function plotAxes(x_lr_norm, y_ap_norm, z_dv_norm)
    %PLOTAXES Plot axis vectors
    % Plot the vectors
    figure;
    x_norm_mat = [x_lr_norm(1); y_ap_norm(1); z_dv_norm(1)];
    y_norm_mat = [x_lr_norm(2); y_ap_norm(2); z_dv_norm(2)];
    z_norm_mat = [x_lr_norm(3); y_ap_norm(3); z_dv_norm(3)];
    pos_mat = [0; 0; 0];
    quiver3(pos_mat, pos_mat, pos_mat, x_norm_mat, y_norm_mat, z_norm_mat)
    text(x_lr_norm(1), x_lr_norm(2), x_lr_norm(3), "u1 (AP)");
    text(y_ap_norm(1), y_ap_norm(2), y_ap_norm(3), "u2");
    text(z_dv_norm(1), z_dv_norm(2), z_dv_norm(3), "u3");
    xlabel('X-axis')
    ylabel('Y-axis')
    zlabel('Z-axis')
end

