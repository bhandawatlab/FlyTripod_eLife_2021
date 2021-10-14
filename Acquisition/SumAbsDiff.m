function sad_value = SumAbsDiff(input_image, ref_image)
    %SUMABSDIFF Summary of this function goes here
    %   Detailed explanation goes here
%     numRows = size(mInputImage, 1);
%     numCols = size(mInputImage, 2);

%     blockLength = size(mRefBlock, 1);
%     blockRadius = (blockLength - 1) / 2;

%     mInputImagePadded = padarray(mInputImage, [blockRadius, blockRadius], 'replicate', 'both');
% 
%     mBlockCol = im2col(mInputImagePadded, [blockLength, blockLength], 'sliding');
    % Convert the images to grayscale
    input_image_gray = im2gray(input_image);
    ref_image_gray = im2gray(ref_image);
    
    % Smooth
    input_image_smooth = imgaussfilt(input_image_gray, 8);
    ref_image_smooth = imgaussfilt(ref_image_gray, 8);
    
    % Do edge detection
    input_image_edges = edge(input_image_smooth, 'Roberts');
    ref_image_edges = edge(ref_image_smooth, 'Roberts');
    figure;
    imshow(input_image_edges)
    
    % Take the sum of absolute differences
    sad_value = sum(imabsdiff(input_image_edges, ref_image_edges), 'all');
    fprintf("SAD value: %d\n", sad_value)
    %mSumAbsDiff = col2im(mSumAbsDiff, [blockLength, blockLength], [(numRows + blockLength - 1), (numCols + blockLength - 1)]);
end
