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

global total_frame_count savecount vid src frame_queue diffarray rROI frame_queue_length previous_frame;
frame_rate = 222.22222;  % Hz (Default in Pylon viewer)
recording_duration = 5; % Duration (s)

%Counting total number of preview frames in real time. We need to keep
%track of number of frames to reset the culmulating data every predefined
%number of frames.
total_frame_count = 1;

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
previous_frame = zeros(rROI(4), rROI(3), 1);

%Need this variable to keep track of total number of videos saved.
savecount = 1;

%% Camera 1 setup
[vid, src] = setDefaultCameraParameters(1, frame_rate, recording_duration);

%% Camera 2 setup
%Creating video and source object
[vid2, src2] = setDefaultCameraParameters(2, frame_rate, recording_duration);

%% Create the preview window
% Create the figure
preview_figure = figure;

%% Set the ROI for camera 1
% Load the ROI
default_camera1_roi_data = load('D:\GitHub\FlyTripod_eLife_2021\Acquisition\Camera1_ROI.mat', 'camera1_roi');
default_roi = default_camera1_roi_data.camera1_roi.Position;
% default_roi = [0,0,1000,1000];

% Start video streaming
subplot(2,2,1);  % Set the subplot section
blank_image_data = zeros(rROI(4), rROI(3), 3, 'uint8');  % 1264x1984 uint8
camera1_image_obj = imshow(blank_image_data);
axis image;
preview(vid, camera1_image_obj);
start(vid);

% Start the ROI selection
camera1_roi = images.roi.Rectangle(gca,'Position', default_roi, 'FaceAlpha', 0);

% Pause and save the ROI
pause;
save('D:\GitHub\FlyTripod_eLife_2021\Acquisition\Camera1_ROI.mat', 'camera1_roi');

% Update and close the camera ROI
stop(vid);
vid.ROIPosition = camera1_roi.Position;
delete(camera1_roi);
start(vid);

%% Set the ROI for camera 2
% Load the ROI
default_camera2_roi_data = load('D:\GitHub\FlyTripod_eLife_2021\Acquisition\Camera2_ROI.mat', 'camera2_roi');
default_roi = default_camera2_roi_data.camera2_roi.Position;

% Start video streaming
subplot(2,2,2);  % Set the subplot section
blank_image_data = zeros(rROI(4), rROI(3), 3, 'uint8');  % 1264x1984 uint8
camera2_image_obj = imshow(blank_image_data);
axis image;
p = preview(vid2, camera2_image_obj);
start(vid2);

% Start the ROI selection
camera2_roi = images.roi.Rectangle(gca,'Position', default_roi, 'FaceAlpha', 0);

% Pause and save the ROI
pause;
save('D:\GitHub\FlyTripod_eLife_2021\Acquisition\Camera2_ROI.mat', 'camera2_roi');

% Update and close the camera ROI
stop(vid2);
vid2.ROIPosition = camera2_roi.Position;
delete(camera2_roi);
start(vid2);

% This allows you to share data between graphics objects. It's a clever way
% of grabbing the image frame whenever the camera updates the preview
% window object 'hImage'. This update is set up below in 'preview(vid,
% hImage)' which states that whenever the camera updates the preview, to
% call the function 'PreviewUpdateD' to do something with this new frame.
% Only observe one camera to avoid conflicts.
setappdata(camera1_image_obj,'UpdatePreviewWindowFcn',@PreviewUpdate);

%%
%-----------------Acquisition Starts-----------------%

% preview(vid, hImage);
% start(vid);
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
