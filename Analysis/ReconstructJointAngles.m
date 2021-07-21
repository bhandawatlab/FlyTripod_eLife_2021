%% Reconstruct joint angles from 3D leg positions in the walking chamber
clear all; close all;

% Frame rate (FPS). This is the resulting frame rate after setting it to
% 387 FPS
frame_rate = 374.95;

% Choose the date and trial number
recording_date = datetime('3/21/2018');
trial_number = 266;

% Input the start and end frames to analyze (where walking occurs and
% tracking is good) [1-indexed: Add 1 if 0-indexed from QuickTime player]
good_frames = [204 453];

% Set the cutoff for the tracking likelihood (DLC default is 0.6). These
% values will be set to NaN and ignored.
pcutoff = 0.6;

% Set the max speed cut-off for ignoring tracking errors
max_speed = 80;  % mm/s

% Read the data file
data_file = "D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Jonathan\Walking Chamber analysis of existing data\DataFiles.xlsx";
T=readtable(data_file);

% Get all data rows with the given parameters
trial_numbers = str2num(char(T.trial_number));
recording_dates = T.recording_date;
fly_data = T(and(recording_dates == recording_date, trial_numbers == trial_number), :);

% Load the pixel spacing
pixel_spacing_data = load("D:\GitHub\FlyTripod_eLife_2021\Preprocessing\PixelSpacing.mat", '-mat');
mm_per_pixel = pixel_spacing_data.mm_pixel;

%% Run 3D reconstruction
% Update the tracking data file with 3D-reconstructed points
tracking_data_file = fly_data.tracking_file{1};
tracking_data = updateFly3DTrackingData(tracking_data_file, mm_per_pixel, frame_rate, max_speed, pcutoff);

%% Plot the XYZ points for each joint
close all
plotJointXYZPoints(tracking_data);

% %% Plot the limb lengths
% limb_length_data = plotLimbLengths(tracking_data);

% %% Save the 3D leg joint positions as a movie
% save3DJointMovie(tracking_data_file, tracking_data);

% % 3D plot
% figure
% plot3(CTr_R_Pro_XYZ_mm(:,1),CTr_R_Pro_XYZ_mm(:,2),CTr_R_Pro_XYZ_mm(:,3))
