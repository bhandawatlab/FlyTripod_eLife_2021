function [cam1_filepath, cam2_filepath] = generateVideoFilepaths()
    %generateVideoFilepath Summary of this function goes here
    %   Detailed explanation goes here
    global output_folder;
    date_str = datestr(date, 'yymmdd');
    date_subfolder = fullfile(output_folder, date_str);
    video_index = 1;
    cam1_filename = sprintf("%s_cam1_%d.avi", date_str, video_index);
    cam1_filepath = fullfile(date_subfolder, cam1_filename);
    cam2_filename = sprintf("%s_cam2_%d.avi", date_str, video_index);
    cam2_filepath = fullfile(date_subfolder, cam2_filename);
    if isfolder(date_subfolder)
        while isfile(cam1_filepath) || isfile(cam2_filepath)
            video_index = video_index + 1;
            cam1_filename = sprintf("%s_cam1_%d.avi", date_str, video_index);
            cam1_filepath = fullfile(date_subfolder, cam1_filename);
            cam2_filename = sprintf("%s_cam2_%d.avi", date_str, video_index);
            cam2_filepath = fullfile(date_subfolder, cam2_filename);
        end
    else
        mkdir(date_subfolder)
    end
end
