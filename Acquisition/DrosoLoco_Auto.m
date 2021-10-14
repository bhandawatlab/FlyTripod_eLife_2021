%This code initiates automatic video acquisition.
%Make sure to calibrate thresholds using the codes in
%"/ThresholdCalibrator".
%
%How this code works:
% This code will call PreviewUpdate.m, which can extract current
% frame data from the camera. Sum of the pixel intensity values of the
% current frame will be taken. This sum will then be compared with that of
% the previous frame, by subtraction. If the output of subtracion is above
% a predefined threshold (variable t1), this code assumes that a motion was
% detected. This initiates 4 seconds of video acquisition (variable record
% time). During this period, the preview window will freeze.
% After the acquisition, the video in memory will be quickly analyzed to
% check if the video is worth being saved to disk. This is judged by
% checking whether the fly walked across a certain amount of space
% (function iswalking.m).
%
%Before running:
% run the DrosoAutoCalibrator.m in /ThresholdCalibrator to see what
% threshold values should be used. After determining good values for the 
% thresholds, assign them to variable t1 in function PreviewUpdate.m and
% variable threshold in function iswalking.m
%
% @Chanwoo Chun <cc2465@cornell.edu>

% Data gets saved here: D:\Documents\walking_project\Data\RawData

clear; close all; daqreset;
imaqreset;

global total_frame_count savecount vid src vid2 src2 frame_queue diffarray rROI frame_queue_length previous_frame previous_sad;
frame_rate = 222.22222;  % Hz (Default in Pylon viewer)
recording_duration = 5; % Duration (s)

%Counting total number of preview frames in real time. We need to keep
%track of number of frames to reset the culmulating data every predefined
%number of frames.
total_frame_count = 1;
full_roi = [0, 0, 1984, 1264];  % Full ROI for reference

% Array that holds 10 frames for use in averaging and determining motion
% whenever a new frame is acquired (This previously held 1000 frames but we
% get the following error: [Error using zeros: Requested 1264x1984x1000
% (18.7GB) array exceeds maximum array size preference. Creation of arrays
% greater than this limit may take a long time and cause MATLAB to become
% unresponsive.)]

% TODO: This frame queue is too large to be manageable with high resolution
% data. Change this to a single value (the motion score) or else we won't
% be able to save data for longer than a couple of seconds.
% frame_queue = zeros(rROI(4), rROI(3), frame_queue_length);
previous_frame = NaN;

%Need this variable to keep track of total number of videos saved.
savecount = 1;

%% Camera 1 setup
[vid, src] = setDefaultCameraParameters(1, frame_rate, recording_duration);

%% Camera 2 setup
%Creating video and source object
[vid2, src2] = setDefaultCameraParameters(2, frame_rate, recording_duration);

%% Create the preview window
% Create the figure
preview_figure = figure('units','normalized','outerposition',[0 0 1 1]);

%% Set the ROI for camera 1
% Load the ROI
default_camera1_roi_data = load('D:\GitHub\FlyTripod_eLife_2021\Acquisition\Camera1_ROI.mat', 'camera1_roi');
camera1_roi = default_camera1_roi_data.camera1_roi;
% default_roi = [0,0,1000,1000];

% Start video streaming
subplot(1,2,1);  % Set the subplot section
blank_image_data = zeros(full_roi(4), full_roi(3), 3, 'uint8');  % 1264x1984 uint8
camera1_image_obj = imshow(blank_image_data);
axis image;
title("Camera 1")
preview(vid, camera1_image_obj);
camera1_axis = gca;
start(vid);

% Start the ROI selection
camera1_roi_data = images.roi.Rectangle(gca,'Position', camera1_roi, 'FaceAlpha', 0);

% Pause and save the ROI
pause;
try
    updated_roi = camera1_roi_data.Position;
    save('D:\GitHub\FlyTripod_eLife_2021\Acquisition\Camera1_ROI.mat', 'camera1_roi');
    
    % Update and close the camera ROI
    stop(vid);
    camera1_roi = updated_roi;
    vid.ROIPosition = camera1_roi;
    delete(camera1_roi_data);
    start(vid);
catch ME
    % Some error occurred if you get here.
    errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message);
    fprintf(1, '%s\n', errorMessage);
    uiwait(warndlg(errorMessage));
end

%% Set the ROI for camera 2
% Load the ROI
default_camera2_roi_data = load('D:\GitHub\FlyTripod_eLife_2021\Acquisition\Camera2_ROI.mat', 'camera2_roi');
camera2_roi = default_camera2_roi_data.camera2_roi;

% Start video streaming
subplot(1,2,2);  % Set the subplot section
blank_image_data = zeros(full_roi(4), full_roi(3), 3, 'uint8');  % 1264x1984 uint8
camera2_image_obj = imshow(blank_image_data);
axis image;
title("Camera 2")
preview(vid2, camera2_image_obj);
camera2_axis = gca;
start(vid2);

% Start the ROI selection
camera2_roi_data = images.roi.Rectangle(gca,'Position', camera2_roi, 'FaceAlpha', 0);

% Pause and save the ROI
pause;
try
    updated_roi = camera2_roi_data.Position;
    save('D:\GitHub\FlyTripod_eLife_2021\Acquisition\Camera2_ROI.mat', 'camera2_roi');
    
    % Update and close the camera ROI
    stop(vid2);
    camera2_roi = updated_roi;
    vid2.ROIPosition = camera2_roi;
    delete(camera2_roi_data);
    start(vid2);
catch ME
    % Some error occurred if you get here.
    errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message);
    fprintf(1, '%s\n', errorMessage);
    uiwait(warndlg(errorMessage));
end

% This allows you to share data between graphics objects. It's a clever way
% of grabbing the image frame whenever the camera updates the preview
% window object 'hImage'. This update is set up below in 'preview(vid,
% hImage)' which states that whenever the camera updates the preview, to
% call the function 'PreviewUpdateD' to do something with this new frame.
% Only observe one camera to avoid conflicts.
setappdata(camera2_image_obj,'UpdatePreviewWindowFcn',@PreviewUpdate);

%%
%-----------------Acquisition Starts-----------------%
% Wait for user input
disp('Press any key to stop.');
pause;

disp('Stopped');

% Stop camera 1
stoppreview(vid);
delete(vid);
clear vid;

% Stop camera 2
stoppreview(vid2);
delete(vid2);
clear vid2;

% Close the figure
delete(preview_figure);
close all; 
daqreset; imaqreset; clear all;
