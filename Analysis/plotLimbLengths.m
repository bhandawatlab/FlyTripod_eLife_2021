function limb_lengths = plotLimbLengths(tracking_data)
    %PLOTLIMBLENGTHS Plot a summary of the limb lengths for a single fly.
    % Set up a structure describing the proximal and distal joint field
    % names for each limb
    femur_field = 'femur_joints';  femur_joints = {'CTr_R_Pro_XYZ', 'FTi_R_Pro_XYZ'};
    tibia_field = 'tibia_joints';  tibia_joints = {'FTi_R_Pro_XYZ', 'TiTa_R_Pro_XYZ'};
    tarsus_field = 'tarsus_joints';  tarsus_joints = {'TiTa_R_Pro_XYZ', 'Ta_XYZ'};
    limb_joint_fieldnames = struct(...
        femur_field, femur_joints,...
        tibia_field, tibia_joints,...
        tarsus_field, tarsus_joints...
        );
    limb_fieldnames = fieldnames(limb_joint_fieldnames);
    limb_count = length(limb_fieldnames);
    
    % Set up the output structure
    limb_lengths = struct();
    
    % Access the good frame indices
    good_frame_indices = find(~tracking_data.bad_frames);
    
    % Find and plot all limb lengths
    figure;
    for limb_index = 1:limb_count
        limb_fieldname = limb_fieldnames{limb_index};
        proximal_fieldname = limb_joint_fieldnames(1).(limb_fieldname);
        distal_fieldname = limb_joint_fieldnames(2).(limb_fieldname);
        
        % Get the proximal joint positions
        proximal_joint_array = tracking_data.(proximal_fieldname);
        
        % Use the good frames only
        proximal_joint_array = proximal_joint_array(good_frame_indices, :);
        
        % Get the distal joint positions
        distal_joint_array = tracking_data.(distal_fieldname);
        
        % Use the good frames only
        distal_joint_array = distal_joint_array(good_frame_indices, :);
        
        % Get the distance between both joints
        x1 = proximal_joint_array(:, 1);
        y1 = proximal_joint_array(:, 2);
        z1 = proximal_joint_array(:, 3);
        x2 = distal_joint_array(:, 1);
        y2 = distal_joint_array(:, 2);
        z2 = distal_joint_array(:, 3);
        joint_distance_mm =...
            sqrt((x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2);
        
        % Store the limb length
        split_field_str = split(limb_fieldname, '_joints');
        dist_fieldname = [split_field_str{1} '_length_mm'];
        limb_lengths(1).(dist_fieldname) = joint_distance_mm;
        
        % Get the limb length mean and standard deviation in um
        limb_length_mean = mean(joint_distance_mm*1000);
        limb_length_std_dev = std(joint_distance_mm*1000);
        
        % Plot the mean and standard deviation
        errorbar(limb_index, limb_length_mean, limb_length_std_dev, 'x')
        hold on
    end
    
    % Finalize the plot
    hold off
    xlim([0 limb_count+1])
    xtickangle(45)
    xticks(1:limb_count)
    xticklabels({'femur', 'tibia', 'tarsus'})
    ylabel('Length (\mum)')
    title('Limb length mean and std.')
end
