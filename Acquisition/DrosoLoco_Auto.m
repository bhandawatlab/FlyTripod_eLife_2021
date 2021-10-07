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

global total_frame_count savecount vid src frame_queue diffarray rROI frame_queue_length;
frame_rate = 50;  % Hz

%Region I want to record.
% rROI = [416 423 960 520]; %[480 450 1024 736];
rROI = [0, 0, 1984, 1264];

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
recording_duration = 5; % Duration (s)
frame_queue_length = 10 * recording_duration;

% TODO: This frame queue is too large to be manageable with high resolution
% data. Change this to a single value (the motion score) or else we won't
% be able to save data for longer than a couple of seconds.
frame_queue = zeros(rROI(4), rROI(3), frame_queue_length);

%Array whose elements are total pixel intensity difference between two
%consecutive frames.
diffarray = zeros(1, frame_queue_length);

%Need this variable to keep track of total number of videos saved.
savecount = 1;

%Creating video and source object
vid = videoinput('gentl', 1, 'Mono8');
set(vid, 'Timeout', 500)

src = getselectedsource(vid);

%%
%-----------------initialization-----------------%

vid.ReturnedColorSpace = 'grayscale';

vid.ROIPosition = rROI; %ROI for recording

% src.ExposureTime = 4000;
src.ExposureTime = 500;  % [JP] Decrease

src.AcquisitionFrameRateEnable = 'True';
src.AcquisitionFrameRate = frame_rate;

%Very important. The below code may not work. If not, 
%src.SensorReadoutMode = Normal;

vid.LoggingMode = 'memory';

triggerconfig(vid,'Manual');

%recordtime = 15; %second
% recordtime = 0.5; %second [JP] Decrease

vid.FramesPerTrigger = src.AcquisitionFrameRate*recording_duration;

src.SensorReadoutTime

%%
%-----------------Ask for confirmation-----------------%
set(vid,'Timeout',recording_duration+5);
get(vid);
get(src);
display(['Resulting Frame Rate = ', num2str(src.ResultingFrameRate)]);
display(['Acquisition Frame Rate = ', num2str(src.AcquisitionFrameRate)]);

disp('Press any key to confirm the setting');
pause;

%Calculate the ideal camera running time
ideal_elaps = vid.FramesPerTrigger/src.AcquisitionFrameRate;

%----------------preview setup----------------------%

% An image object of the same size as the video is used to store and
% display incoming frames.

% Retrieve the video resolution. 
%vidRes = vid.VideoResolution;
vidRes = vid.ROIPosition;
% Create a figure and an image object.
f=figure('Visible', 'off');

% The Video Resolution property returns values as width by height, but
% MATLAB images are height by width, so flip the values.
% imageRes = fliplr([rROI(3) 1264/2-rROI(2)]);

%The top subplot shows actual live video.
subplot(2,1,1);
blank_image_data = zeros(rROI(4), rROI(3), 3, 'uint8');  % 1264x1984 uint8
hImage = imshow(blank_image_data);
% Set the axis of the displayed image to maintain the aspect ratio of the 
% incoming frame.
axis image;

% This allows you to share data between graphics objects. It's a clever way
% of grabbing the image frame whenever the camera updates the preview
% window object 'hImage'. This update is set up below in 'preview(vid,
% hImage)' which states that whenever the camera updates the preview, to
% call the function 'PreviewUpdateD' to do something with this new frame.
setappdata(hImage,'UpdatePreviewWindowFcn',@PreviewUpdate);

%%
%-----------------Acquisition Starts-----------------%

preview(vid, hImage);
start(vid);
disp('Press any key to stop.');


pause;

disp('Stopped');

stoppreview(vid);
delete(f);

delete(vid);
clear vid;

close all; 
daqreset; imaqreset; clear all;
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

global total_frame_count savecount vid src frame_queue diffarray rROI frame_queue_length;
frame_rate = 50;  % Hz

%Region I want to record.
% rROI = [416 423 960 520]; %[480 450 1024 736];
rROI = [0, 0, 1984, 1264];

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
recording_duration = 5; % Duration (s)
frame_queue_length = 10 * recording_duration;

% TODO: This frame queue is too large to be manageable with high resolution
% data. Change this to a single value (the motion score) or else we won't
% be able to save data for longer than a couple of seconds.
frame_queue = zeros(rROI(4), rROI(3), frame_queue_length);

%Array whose elements are total pixel intensity difference between two
%consecutive frames.
diffarray = zeros(1, frame_queue_length);

%Need this variable to keep track of total number of videos saved.
savecount = 1;

%Creating video and source object
vid = videoinput('gentl', 1, 'Mono8');
set(vid, 'Timeout', 500)

src = getselectedsource(vid);

%%
%-----------------initialization-----------------%

vid.ReturnedColorSpace = 'grayscale';

vid.ROIPosition = rROI; %ROI for recording

% src.ExposureTime = 4000;
src.ExposureTime = 500;  % [JP] Decrease

src.AcquisitionFrameRateEnable = 'True';
src.AcquisitionFrameRate = frame_rate;

%Very important. The below code may not work. If not, 
%src.SensorReadoutMode = Normal;

vid.LoggingMode = 'memory';

triggerconfig(vid,'Manual');

%recordtime = 15; %second
% recordtime = 0.5; %second [JP] Decrease

vid.FramesPerTrigger = src.AcquisitionFrameRate*recording_duration;

src.SensorReadoutTime

%%
%-----------------Ask for confirmation-----------------%
set(vid,'Timeout',recording_duration+5);
get(vid);
get(src);
display(['Resulting Frame Rate = ', num2str(src.ResultingFrameRate)]);
display(['Acquisition Frame Rate = ', num2str(src.AcquisitionFrameRate)]);

disp('Press any key to confirm the setting');
pause;

%Calculate the ideal camera running time
ideal_elaps = vid.FramesPerTrigger/src.AcquisitionFrameRate;

%----------------preview setup----------------------%

% An image object of the same size as the video is used to store and
% display incoming frames.

% Retrieve the video resolution. 
%vidRes = vid.VideoResolution;
vidRes = vid.ROIPosition;
% Create a figure and an image object.
f=figure('Visible', 'off');

% The Video Resolution property returns values as width by height, but
% MATLAB images are height by width, so flip the values.
% imageRes = fliplr([rROI(3) 1264/2-rROI(2)]);

%The top subplot shows actual live video.
subplot(2,1,1);
blank_image_data = zeros(rROI(4), rROI(3), 3, 'uint8');  % 1264x1984 uint8
hImage = imshow(blank_image_data);
% Set the axis of the displayed image to maintain the aspect ratio of the 
% incoming frame.
axis image;

% This allows you to share data between graphics objects. It's a clever way
% of grabbing the image frame whenever the camera updates the preview
% window object 'hImage'. This update is set up below in 'preview(vid,
% hImage)' which states that whenever the camera updates the preview, to
% call the function 'PreviewUpdateD' to do something with this new frame.
setappdata(hImage,'UpdatePreviewWindowFcn',@PreviewUpdate);

%%
%-----------------Acquisition Starts-----------------%

preview(vid, hImage);
start(vid);
disp('Press any key to stop.');


pause;

disp('Stopped');

stoppreview(vid);
delete(f);

delete(vid);
clear vid;

close all; 
daqreset; imaqreset; clear all;
