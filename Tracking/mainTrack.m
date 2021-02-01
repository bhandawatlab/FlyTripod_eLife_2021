%This code starts automatic video processing. Once the video is processed,
%shot.mat file will be saved. This file includes CoM kinematics data and
%leg positions, etc.
%
%@Chanwoo Chun <cc2465@cornell.edu>

close all

matdir = pwd;

%Define the location of the video
rootdir = ['..' filesep '..' filesep 'Data' filesep 'RawVideoExample' filesep '171231Video'];
% variable info is an array whose elements are:
% info(1): video number.
% info(2): first frame to begin analysis.
% info(3): final frame to finish analysis.
info = [165	500	700];
dataNum = info(1);
frBegin = info(2);
frEnd = info(3);

%Grab the video
M=load([rootdir filesep num2str(dataNum)]);

%store the first&last frame number
shot.dataNum = dataNum;
shot.timeStamp = M.data.timeStamp(frBegin:frEnd);
shot.frame = cumsum([frBegin; ones((frEnd-frBegin),1)]);

%get background
background = trueBackground(squeeze(M.data.frames));

%get initial CoM position
shot.com = getCOM3(M,frBegin,frEnd,background);

%binarize the bottom view portion of the video.
[binVid, legTips, notTouching, addedFirst, addedLast] = bottomViewAnalysis(M,frBegin,frEnd,shot.com.worldinfo.wallUpper,shot.com.worldinfo.wallLower,background);
shot.validity.notTouching=notTouching;

%determine validity of the frames. Here, by valid, we mean 1) no wall
%touching, 2) small slip angle (beta), 3) small turn angle (yaw).
[smallBeta, smallYaw, valid] = validity(shot);

shot.validity.smallBeta=smallBeta;
shot.validity.smallYaw=smallYaw;
shot.validity.valid=valid;


[~,name,~] = fileparts(rootdir);
analysisdir = ['..' filesep '..' filesep 'Data' filesep 'ProcessedData' filesep name];
mkdir(analysisdir);

%leg tracking
shot.leg = getLegs2(shot, binVid, legTips, analysisdir, addedFirst, addedLast);

%save processed data
save([analysisdir filesep 'shot_EXAMPLE.mat'],'shot');
