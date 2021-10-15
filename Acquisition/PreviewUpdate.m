% This callback function updates the displayed frame and perform analysis
% on current frame.
% @Chanwoo Chun, <cc2465@cornell.edu>

function PreviewUpdate(obj,event,hImage)

global total_frame_count savecount vid vid2 previous_frame previous_sad trigger_obj daq_output sample_count;

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
        if length(previous_sad) == 100
            sad_mean = mean(previous_sad);
            if sad_mean > sad_threshold
                fprintf("Motion detected: %.1f\n", sad_mean)
                disp('Recording...');
                
                %% Start acquisition
                % Stop the videos to update the trigger settings
                stop(vid);
                stop(vid2);
                
                % Set the trigger modes to external
                src.TriggerMode = 'on';
                src2.TriggerMode = 'on';
                triggerconfig(vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
                triggerconfig(vid2, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
                
                % Update the source trigger parameters
                vid.FramesPerTrigger = 1;
                vid.TriggerRepeat = sample_count - 1;
                src1 = getselectedsource(vid);
                src1.TriggerSelector = 'FrameStart';
                src1.TriggerSource = 'Line1';
                src1.TriggerActivation = 'RisingEdge';
                src1.TriggerMode = 'on';
                vid2.FramesPerTrigger = 1;
                vid2.TriggerRepeat = sample_count - 1;
                src2 = getselectedsource(vid2);
                src2.TriggerSelector = 'FrameStart';
                src2.TriggerSource = 'Line1';
                src2.TriggerActivation = 'RisingEdge';
                src2.TriggerMode = 'on';
                
                % Start the videos again
                start(vid);
                start(vid2);
                
                %trigger(vid);
%                 trigger_obj.preload(daq_output);
                write(trigger_obj, daq_output);  % https://www.mathworks.com/help/daq/daq.interfaces.dataacquisition.write.html
                disp('Recording complete.');
                
                %% Save the camera 1 video
                disp('Saving camera 1...');
                [frames, timeStamp] = getdata(vid);
                s = struct();
                s(1).Camera1Timestamps = timeStamp;
                s(1).Camera1Frames = frames;
                % savecount=savecount+1;
                clear frames;
                
                %% Save the camera 2 video
                disp('Saving camera 2...');
                [frames, timeStamp] = getdata(vid2);
                s(1).Camera2Timestamps = timeStamp;
                s(1).Camera2Frames = frames;
                saver(s);
                savecount=savecount+1;
                clear frames;
                %start(vid2);
                disp('Saving complete.');
                
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
