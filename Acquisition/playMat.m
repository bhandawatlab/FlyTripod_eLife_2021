%This code converts the acquired .mat video files into compressed .avi files.
%@Chanwoo Chun, <cc2465@cornell.edu>

%defined below is a example folder location
rootdir = ['..' filesep '..' filesep 'Data' filesep 'ProcessedData' filesep '181213Sideview'];

start=56;    %enter the FIRST file to convert
endd =131;   %enter the LAST file to convert

cd(rootdir);
for filen=start:endd
    M=load(num2str(filen));
   
    diskLogger = VideoWriter(num2str(filen), 'Motion JPEG AVI');
    open(diskLogger);
    numFrames = size(M.data.frames, 4);
    for ii = 1:numFrames
        writeVideo(diskLogger, M.data.frames(:,:,:,ii));
    end
    close(diskLogger);
    
end

clear all;