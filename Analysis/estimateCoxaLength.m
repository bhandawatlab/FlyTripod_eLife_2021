function coxa_length = estimateCoxaLength(limb_lengths)
    %ESTIMATECOXALENGTH Estimate the coxa length.
    % Use the proportion among coxa-femur-tibia segment lengths to find the
    % coxa lengths. These proportions are from Arena et al., 2012
    % (https://ieeexplore.ieee.org/abstract/document/6290809)
    % Front leg: [0.1, 0.3793, 0.46]
    % Middle leg: [0.1, 0.4167, 0.6]
    % Back leg: [0.1, 0.4065, 0.67]
    
    % We are only looking at the front (prothoracic) leg for now:
    coxa_femur_tibia_proportion = [0.2174 0.8246 1.0000];  % Dividing by the max
    coxa_proportion_of_tibia = coxa_femur_tibia_proportion(1);
    
    % Get the mean tibia length
    tibia_length = limb_lengths.('tibia_mean_length_um');

    % Estimate the coxa length (um)
    coxa_length = tibia_length * coxa_proportion_of_tibia;
end
