%This code will get the length of the middle leg
clear all
close all

rootdir = 'E:\CHUN\ImageCollect\Data\180924Sideview';
filen = 19;
M=load([rootdir filesep num2str(filen)]);

data = M.data.frames*2;

%start with a tarsi (foot) of the middle leg in the bottom view.
%After done with registering joints in the bottom view, move on to the side
%view. Press 1 for registration of each joint.
%Leg length will be calculated through 5 frames.
Rn=getLegLengthEditor(data)

close all