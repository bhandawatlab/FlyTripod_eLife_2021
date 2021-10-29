%% Reconstruct joint angles from 3D leg positions only during swing phase
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

%% Plot the leg relative to the AP-axis with the origin at the posterior point
ap_axis_points = getAPAxisPoints(tracking_data);
% plotAPAxisInterpolatedLeg(ap_axis_points, tracking_data_file);

%% Get the swing phases from the tarsus speed minima
close all;

% Get the tarsus height (Dorsal-ventral Z-axis)
tarsus_height = ap_axis_points.Ta_U(:, 3);
figure; plot(tarsus_height)
title("Tarsus Height")
hold on

% Plot the baseline height
polynomialOrder = 1;
windowWidth = 69;
baseline_height = sgolayfilt(tarsus_height, polynomialOrder, windowWidth);
% detrendedY = y - baselineY;
plot(baseline_height)
hold off

% Set to NaN all values below the baseline
tarsus_height_swing_phase = tarsus_height;
tarsus_height_swing_phase(tarsus_height <= baseline_height) = NaN;
figure; plot(tarsus_height_swing_phase)
title("Tarsus Height - Baseline")

% Split up each swing phase and store frames in a cell array if it
% satisfies the duration condition
min_swing_duration = 10;  % 10 ms
min_swing_frame_count = (min_swing_duration / 1000) * frame_rate;
swing_phase_frames = {};
current_swing_height = [];
current_swing_frames = [];
for n=1:length(tarsus_height_swing_phase)
    height_value = tarsus_height_swing_phase(n);
    
    % If we have reached the end of a swing, store the current swing
    if isnan(height_value)
        if ~isempty(current_swing_frames)
            % Check if the swing satisfies the duration condition
            if length(current_swing_frames) >= min_swing_frame_count
                
                % Check if the swing satisfies the condition of only having
                % one local maxima
                smoothed_swing_height = smoothdata(current_swing_height);
                %figure; plot(smoothed_swing)
                if length(findpeaks(smoothed_swing_height)) == 1
                    swing_phase_frames{end+1} = current_swing_frames;
                end
            end
        end
        current_swing_frames = [];
        current_swing_height = [];
    elseif ~isnan(height_value)
        current_swing_frames = [current_swing_frames; n];
        current_swing_height = [current_swing_height; height_value];
    end
end

%% Add the swing phase frames to the tracking data file
tracking_data(1).swing_phase_frames = swing_phase_frames;
save(tracking_data_file, '-struct', 'tracking_data');