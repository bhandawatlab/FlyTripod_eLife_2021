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

clear all; close all; daqreset; 
imaqreset;

global prevcount savecount vid src tempframes diffarray rROI;

%Region I want to record.
rROI = [480 450 1024 736];

%Counting total number of preview frames in real time. We need to keep
%track of number of frames to reset the culmulating data every predefined
%number of frames.
prevcount = 1;

%tempframes is declared to be used in other function. This temporarily
%holds current frame image.
tempframes = zeros(1264/2-rROI(2),rROI(3),1000);

%Array whose elements are total pixel intensity difference between two
%consecutive frames.
diffarray = zeros(1, 1000);

%Need this variable to keep track of total number of videos saved.
savecount = 1;

%Creating video and source object
vid = videoinput('gentl', 1, 'Mono8');
src = getselectedsource(vid);

%%
%-----------------initialization-----------------%

vid.ReturnedColorSpace = 'grayscale';

vid.ROIPosition = rROI; %ROI for recording

src.ExposureTime = 2500;

src.AcquisitionFrameRate = 387;

%Very important. The below code may not work. If not, 
%src.SensorReadoutMode = Normal;

vid.LoggingMode = 'memory';

triggerconfig(vid,'Manual');

recordtime = 4; %second

vid.FramesPerTrigger = src.AcquisitionFrameRate*recordtime;

src.SensorReadoutTime

%%
%-----------------Ask for confirmation-----------------%
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
imageRes = fliplr([rROI(3) 1264/2-rROI(2)]);

%The top subplot shows actual live video.
subplot(5,1,1);
hImage = imshow(zeros(imageRes));
% Set the axis of the displayed image to maintain the aspect ratio of the 
% incoming frame.
axis image;

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
