function [video_obj, video_src] = setDefaultCameraParameters(camera_index, frame_rate, recording_duration)
    %SETDEFAULTCAMERAPARAMETERS Summary of this function goes here
    %   Detailed explanation goes here%Creating video and source object
    video_obj = videoinput('gentl', camera_index, 'Mono8');
    set(video_obj, 'Timeout', 500)
    video_src = getselectedsource(video_obj);
    video_obj.ReturnedColorSpace = 'grayscale';
    % vid.ROIPosition = rROI; %ROI for recording
    % src.ExposureTime = 4000;
    video_src.ExposureTime = 500;  % [JP] Decrease
    video_src.AcquisitionFrameRateEnable = 'True';
    video_src.AcquisitionFrameRate = frame_rate;

    %Very important. The below code may not work. If not, 
    %src.SensorReadoutMode = Normal;
    video_obj.LoggingMode = 'memory';
    triggerconfig(video_obj,'Manual');
    video_obj.FramesPerTrigger = video_src.AcquisitionFrameRate*recording_duration;
    video_src.SensorReadoutTime
    set(video_obj,'Timeout',recording_duration+5);
    
    % Show the camera parameters
    % get(vid);
    % get(src);
    display(['(Camera ' num2str(camera_index) ') Resulting Frame Rate = ', num2str(video_src.ResultingFrameRate)]);
end

