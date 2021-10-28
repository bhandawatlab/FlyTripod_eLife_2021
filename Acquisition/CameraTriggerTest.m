output_file = 'D:\GitHub\FlyTripod_eLife_2021\Acquisition\TestVideo2.avi';
frameRate = 100; %20; %  20; %200;  %320;
numSecs = 10; % Ephys. recording duration
numFrames = ceil(numSecs * frameRate);

% Grab the camera
vid1 = videoinput('gentl', 1, 'Mono8');  % For IR lights
config = triggerinfo(vid1);
triggerconfig(vid1, 'hardware', 'DeviceSpecific', 'DeviceSpecific'); %configure trigger type that cameras will expect

% Set up the video trigger
vid1.FramesPerTrigger = 1;
vid1.TriggerRepeat = numFrames - 1;
vid1.LoggingMode = 'Disk';
vid1_writer = VideoWriter(output_file, 'Uncompressed AVI');
vid1_writer.FrameRate = frameRate;
vid1.DiskLogger = vid1_writer;
src1 = getselectedsource(vid1);
src1.TriggerSelector = 'FrameStart';
src1.TriggerSource = 'Line1';
src1.TriggerActivation = 'RisingEdge';
src1.TriggerMode = 'on';
src1.AcquisitionFrameRateEnable = 'False';
src1.AcquisitionFrameRate = frameRate;
src1.DeviceLinkThroughputLimitMode = 'Off';

% Set up the DAQ output
s = daq('ni');
s.Rate = frameRate*2;
s.addoutput('Dev3','ao0','voltage'); % Dual Basler camera trigger    

%% Start the normal ephys. acquisition
disp('Data acquisition starting.');

% Set up the DAQ output signal
ephys_output = zeros(numFrames*2, 1);
ephys_output(1:2:end) = 10;

% Start the camera
disp('Starting the camera');
start(vid1);

% Load the voltage trigger
s.preload(ephys_output);
disp(['Number queued: ' num2str(s.NumScansQueued)]);

% Start the acquisition
disp('Starting acquisition')
s.start();
while s.Running
    pause(0.5)
    fprintf("While loop: Frames acquired = %d\n", vid1.FramesAcquired)
end

fprintf("Acquisition completed with %d frames acquired\n", vid1.FramesAcquired);
disp(['Frames acquired (vid1): ' num2str(vid1.FramesAcquired)]);
disp(['Frames logged: ' num2str(vid1.DiskLoggerFrameCount)]);
closepreview(vid1);

% Remove image acquisition objects from memory
delete(vid1)
clear vid1

daqreset;
pause(1);
