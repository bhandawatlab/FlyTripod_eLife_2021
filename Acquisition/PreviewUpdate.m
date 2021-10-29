% This callback function updates the displayed frame and perform analysis
% on current frame.
% @Chanwoo Chun, <cc2465@cornell.edu>

function PreviewUpdate(obj,event,hImage)

global frame_rate total_frame_count savecount vid vid2 previous_frame previous_sad trigger_obj daq_output sample_count;

    try
        % Threshold for the sum of absolute differences
        sad_threshold = 5500;

        % Get the frame
        current_frame = event.Data;

        % Display the current image frame.
        set(hImage, 'CData', current_frame);
        
        % Quantify motion in a given range of frames
        previous_motion_score = SumAbsDiff(current_frame, previous_frame);
        previous_sad = [previous_sad previous_motion_score];            
        if length(previous_sad) == 1000
            sad_mean = mean(previous_sad);
            if sad_mean > sad_threshold
                fprintf("Motion detected: %.1f\n", sad_mean)
                disp('Recording...');
                
                %% Start acquisition
                % Stop the videos to update the trigger settings
                stop(vid);
                stop(vid2);
                
                % Generate the video filepaths
                [cam1_filepath, cam2_filepath] = generateVideoFilepaths();
                
                % Set the disk logging parameters
                vid.LoggingMode = 'Disk';
                vid_writer = VideoWriter(cam1_filepath, 'Uncompressed AVI');
                vid_writer.FrameRate = frame_rate;
                vid.DiskLogger = vid_writer;
                
                vid2.LoggingMode = 'Disk';
                vid2_writer = VideoWriter(cam2_filepath, 'Uncompressed AVI');
                vid2_writer.FrameRate = frame_rate;
                vid2.DiskLogger = vid2_writer;
                
                % Set the trigger modes to external
                triggerconfig(vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
                src1 = getselectedsource(vid);
                src1.TriggerMode = 'on';
                vid.FramesPerTrigger = 1;
                vid.TriggerRepeat = sample_count - 1;
                src1.TriggerSelector = 'FrameStart';
                src1.TriggerSource = 'Line1';
                src1.TriggerActivation = 'RisingEdge';
                src1.AcquisitionFrameRateEnable = 'False';
                src1.AcquisitionFrameRate = frame_rate;
                src1.DeviceLinkThroughputLimitMode = 'Off';
                
                triggerconfig(vid2, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
                src2 = getselectedsource(vid2);
                src2.TriggerMode = 'on';
                vid2.FramesPerTrigger = 1;
                vid2.TriggerRepeat = sample_count - 1;
                src2.TriggerSelector = 'FrameStart';
                src2.TriggerSource = 'Line1';
                src2.TriggerActivation = 'RisingEdge';
                src2.AcquisitionFrameRateEnable = 'False';
                src2.AcquisitionFrameRate = frame_rate;
                src2.DeviceLinkThroughputLimitMode = 'Off';
                
                % Start the videos again
                disp('Starting cameras');
                start(vid);
                start(vid2);
                
                % Load the voltage trigger
                trigger_obj.preload(daq_output);
                disp(['Number queued: ' num2str(trigger_obj.NumScansQueued)]);

                % Start the acquisition
                trigger_obj.start();
                while trigger_obj.Running
                    pause(0.5)
                    fprintf("While loop: Frames acquired = %d, %d\n", vid.FramesAcquired, vid2.FramesAcquired)
                end

                fprintf("Acquisition completed with %d, %d frames acquired\n", vid.FramesAcquired, vid2.FramesAcquired)
                fprintf("Frames logged: %d, %d\n", vid.DiskLoggerFrameCount, vid2.DiskLoggerFrameCount)
                fprintf("Cam1: %s\nCam2: %s\n", cam1_filepath, cam2_filepath);
                
                % Stop the videos to update the trigger settings
                stop(vid);
                stop(vid2);
                
                % Set the trigger modes back to manual (allows previewing)
                triggerconfig(vid, 'Manual');
                triggerconfig(vid2, 'Manual');
                
                % Start the videos again
                start(vid);
                start(vid2);
            else
                fprintf("%d\n", sad_mean)
            end
            previous_sad = [];
        end

        % Update the previous frame and frame count
        previous_frame = current_frame;
        total_frame_count=total_frame_count+1;

        if savecount>=600
            stoppreview(vid);
            stoppreview(vid2);
        end
    catch ME
        disp("ruh roh")
        
        % Some error occurred if you get here.
        errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
            ME.stack(1).name, ME.stack(1).line, ME.message);
        fprintf(1, '%s\n', errorMessage);
        uiwait(warndlg(errorMessage));
    end
end
