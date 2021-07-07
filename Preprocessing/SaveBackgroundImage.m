clear all; close all;

% Load a video that has background images (no visible fly) that we can save
% for use for background subtraction.
video_file = "D:\Bhandawatlab_Drexel Dropbox\Bhandawat_Lab_Transfer\Chun\ImageCollect\Data\180325Sideview\74.avi";
videoSource = VideoReader(video_file);
videoPlayer = vision.VideoPlayer();
noise_level = 3;
previous_frame = NaN;

% Allows you to exit with ESC key:
H = uicontrol('Style', 'PushButton', ...
                    'String', 'Break', ...
                    'Callback', 'delete(gcbf)');

frame_to_grab = 78;
frame_index = 0;
while hasFrame(videoSource) && ishandle(H) && frame_index < frame_to_grab
    frame_index = frame_index + 1;
    frame = readFrame(videoSource);
    videoPlayer(frame);
    pause(0.1);
end

release(videoPlayer);

% Save the frame
imwrite(frame, 'D:\GitHub\FlyTripod_eLife_2021\Preprocessing\180325SideView_BG.jpg');
