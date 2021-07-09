%% Reconstruct joint angles from 3D leg positions in the walking chamber
clear all; close all;
frame_rate = 387; % FPS

% Choose the date and trial number
recording_date = datetime('3/21/2018');
trial_number = 266;

% Input the start and end frames to analyze (where walking occurs and
% tracking is good) [1-indexed: Add 1 if 0-indexed from QuickTime player]
good_frames = [204 453];

% Set the cutoff for the tracking likelihood (DLC default is 0.6). These
% values will be set to NaN and ignored.
pcutoff = 0.6;

% Read the data file
data_file = "D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Jonathan\Walking Chamber analysis of existing data\DataFiles.xlsx";
T=readtable(data_file);

% Get all data rows with the given parameters
trial_numbers = str2num(char(T.trial_number));
recording_dates = T.recording_date;
fly_data = T(and(recording_dates == recording_date, trial_numbers == trial_number), :);

% Load the tracking data
tracking_data = load(fly_data.tracking_file{1}, '-mat');

% Load the pixel spacing
pixel_spacing_data = load("D:\GitHub\FlyTripod_eLife_2021\Preprocessing\PixelSpacing.mat", '-mat');
mm_per_pixel = 1/pixel_spacing_data.pixels_mm;

%% Run 3D reconstruction
%% CTr R Pro:
% Top:
CTr_R_Pro_Top_xy = [tracking_data.CTr_R_Pro_Top.x' tracking_data.CTr_R_Pro_Top.y'];

% Ignore all tracks with low likelihood
CTr_R_Pro_Top_lh = tracking_data.CTr_R_Pro_Top.likelihood;
CTr_R_Pro_Top_xy(CTr_R_Pro_Top_lh <= pcutoff, 1:2) = NaN;

% Bottom:
CTr_R_Pro_Bottom_xy = [tracking_data.CTr_R_Pro_Bottom.x' tracking_data.CTr_R_Pro_Bottom.y'];

% Ignore all tracks with low likelihood
CTr_R_Pro_Bottom_lh = tracking_data.CTr_R_Pro_Bottom.likelihood;
CTr_R_Pro_Bottom_xy(CTr_R_Pro_Bottom_lh <= pcutoff, 1:2) = NaN;

% Reconstruct
CTr_R_Pro_XYZ_mm = Mirror3DReconstruction(CTr_R_Pro_Top_xy, CTr_R_Pro_Bottom_xy, mm_per_pixel);

% Ignore all values with speeds > 30mm/s (tracking errors)
CTr_R_Pro_XYZ_speed = getSpeedXYZ(CTr_R_Pro_XYZ_mm, frame_rate);
