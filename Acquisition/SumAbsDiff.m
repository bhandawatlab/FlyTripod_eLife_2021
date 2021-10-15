function sad_value = SumAbsDiff(input_image, ref_image)
    %SUMABSDIFF Summary of this function goes here
    % Convert the images to grayscale
    input_image_gray = im2gray(input_image);
    ref_image_gray = im2gray(ref_image);
    
    % Do edge detection
    input_image_edges = edge(input_image_gray, 'Roberts');
    ref_image_edges = edge(ref_image_gray, 'Roberts');
    
    % Take the sum of absolute differences
    if isnan(ref_image)
        % Skip the first NaN frame
        sad_value = 0;
    else
        sad_value = sum(imabsdiff(input_image_edges, ref_image_edges), 'all');
    end
end
