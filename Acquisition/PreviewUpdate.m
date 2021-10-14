% This callback function updates the displayed frame and perform analysis
% on current frame.
% @Chanwoo Chun, <cc2465@cornell.edu>

function PreviewUpdate(obj,event,hImage)

global total_frame_count savecount frame_queue diffarray vid rROI w frame_queue_length previous_frame;

    try
        % Threshold for the sum of absolute differences
        sad_threshold = 20;

        numavg = 2; %Number of frames to average. Need this to smooth out noise.
        %Change the value to appropriate value for your application, if needed. 

        % Get the frame
        current_frame = event.Data;

        % Display the current image frame.
        set(hImage, 'CData', current_frame);

        % Add it to the frame queue
        current_frame_index = total_frame_count;
%         frame_queue(:,:,current_frame_index)= current_frame;

        % If the number of frames we want to average is reached
        if rem(current_frame_index,numavg)==0
            %Quantify motion in a given range of frames
%             diffarray(current_frame_index-(numavg-1):current_frame_index)=diffavg(frame_queue(:,:,current_frame_index-(numavg-1):current_frame_index));
%             diffarray = diffarray(1:current_frame_index);
%             
%             previous_max_diff = sum(diff(norm(frame_queue(:,:,current_frame_index-(numavg-1):current_frame_index), 1, 3)), 'all');
%             disp("Previous max: ");
%             disp(previous_max_diff);
% 
%             %If the quantified value is above a certain threshold:
%             previous_motion_score = diffarray(current_frame_index);
            previous_motion_score = SumAbsDiff(current_frame, previous_frame);
%             if previous_motion_score >= t1
% 
%                 %start acquisition
%                 trigger(vid);
%                 disp('Recording...');
%                 [frames, timeStamp] = getdata(vid);
%                 s = struct('timeStamp', timeStamp, 'frames', frames);
%                 framesq = squeeze(frames);
% 
%                 %Check if the fly walked across the frame. If so, save the vid.
% %                 if iswalkingD(framesq) == true
%                 current_max_diff = sum(diff(norm(framesq), 1, 3), 'all');
%                 disp("Current max: ");
%                 disp(current_max_diff);
%             
%                 post_acquisition_motion_score = diffavg(framesq);  % Quantify motion in a given range of frames
%                 if post_acquisition_motion_score >= t1
%                     saver(s);
%                     savecount=savecount+1;
%                 else
%                     disp('Not worthy. Resume monitering... Press any key to stop.');
%                     clear frames;
%                 end
%                 start(vid);
%             end
        end

        subplot(2,1,2);
%         plot(diffarray);
        drawnow


        if total_frame_count==frame_queue_length
            total_frame_count=0;
%             clear frame_queue
%             clear diffarray
        end

        % Update the previous frame and frame count
        previous_frame = current_frame;
        total_frame_count=total_frame_count+1;

        if savecount>=600
            stoppreview(vid);
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
