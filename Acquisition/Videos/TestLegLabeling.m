% ventral_video_file = "D:\GitHub\FlyTripod_eLife_2021\Acquisition\Videos\bottom-video.avi";
% frame_range = [26 34];
% labelLeg(ventral_video_file, frame_range);

side_video_file = "D:\GitHub\FlyTripod_eLife_2021\Acquisition\Videos\side-video.avi";
frame_range = [52 56];
labelLeg(side_video_file, frame_range);

disp("Done.")

function labelLeg(video_file, frame_range)
    % Get the filename
    [filepath,name,ext] = fileparts(video_file);

    % Label the leg in video
    v = VideoReader(video_file);
    frames = read(v, frame_range);
    for n=1:abs(frame_range(2) - frame_range(1))
        frame = frames(:,:,:,n);
        imshow(frame)
        
        % Coxa
        hline = gline; % Connect circles
        
        % Femur
        hline1 = gline; % Connect circles
        
        % Tibia
        hline2 = gline; % Connect circles
        uiwait(msgbox('Click when editing is done.'));
        set(hline,'Color','r')
        set(hline1,'Color','r')
        set(hline2,'Color','r')
        
        % Save the frame
        frame_filepath = fullfile(filepath, sprintf("%s_%d.png", name, n));
        saveas(gcf, frame_filepath);
    end
end
