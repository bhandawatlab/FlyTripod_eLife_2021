%Use this code to determine good threshold values for automatic acquisition
%@Chanwoo Chun, <cc2465@cornell.edu>

clear all; close all; daqreset; 
imaqreset;

global prevcount fincount vid B C rROI;

rROI = [480 450 1024 736];%[480 460 1024  779]; %[480 500 1024  615];
prevcount = 1;
B = zeros(1264/2-rROI(2),rROI(3),1000);
C = zeros(1, 1);
fincount = 1;

vid = videoinput('gentl', 1, 'Mono8');
src = getselectedsource(vid);

%%
%-----------------initialization-----------------%

vid.ReturnedColorSpace = 'grayscale';

%rROI = [480 450 1024 736]; %Xoffset, Yoffset, Width, Height

vid.ROIPosition = rROI; %ROI for recording

src.ExposureTime = 2500;

src.AcquisitionFrameRate = 395;


%----------------preview setup----------------------%

% An image object of the same size as the video is used to store and
% display incoming frames.

% Retrieve the video resolution. 
vidRes = vid.ROIPosition;
% Create a figure and an image object.
f=figure('Visible', 'off');

% The Video Resolution property returns values as width by height, but
% MATLAB images are height by width, so flip the values.
imageRes = fliplr([rROI(3) 1264/2-rROI(2)]);

subplot(5,1,1);
hImage = imshow(zeros(imageRes));
%hImage = imshow(zeros(182,1024));
% Set the axis of the displayed image to maintain the aspect ratio of the 
% incoming frame.
axis image;


setappdata(hImage,'UpdatePreviewWindowFcn',@PreviewUpdateC);
%dbtype('PreviewUpdate.m')

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
