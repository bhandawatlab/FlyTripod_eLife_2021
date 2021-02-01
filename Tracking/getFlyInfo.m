%This is a library where a parent function can look up an information on a
%fly data stored in a specific file.
%variable analysisFolder is a path to the file location
%
%@Chanwoo Chun, <cc2465@cornell.edu>

function [gender, strain, weight, legLength, flynum] = getFlyInfo(analysisFolder)

analysisFolder = strrep(analysisFolder,'\',filesep);
[filepath,name] = fileparts(analysisFolder);
vidNum = str2double(name);
[~,name] = fileparts(filepath);
date = str2double(erase(name,'Sideview'));

%% Fly data info
%[info, files] = organize(gender, strain, leglength, expDate, dataRange) %    leg length (mm)  std of leg length (mm) (number of valid steps)
i=1;    
[info{i}, files{i}] = organize(1,'female', 'Berlin K', 2.1032, 180130, 1:4); %     2.1032    0.0898 (3)
i=2; %NO A
[info{i}, files{i}] = organize(2,'female', 'Berlin K', 2.1509, 180130, 9:14); %    2.1509    0.0727 (9)
i=3;
[info{i}, files{i}] = organize(3,'female', 'Berlin K', 2.1840, 180130, 19:90);% >5 steps   2.1840    0.0693 (57)
i=4; %3 and 4 flynum could be a same fly.
[info{i}, files{i}] = organize(4,'female', 'Berlin K', 2.1031, 180130, 109:271);% >5 steps  2.1031    0.0335 (125) 
i=5;
[info{i}, files{i}] = organize(4,'female', 'Berlin K', 2.1031, 180131, 1:144);%Not used for analysis 1.9916    0.0800 (1)
i=6;
[info{i}, files{i}] = organize(5,'female', 'Berlin K', 1.9988, 180131, 147:251); % >5 steps  1.9988    0.0893 (153)
i=7;
[info{i}, files{i}] = organize(6,'female', 'Berlin K', 1.9869, 180131, 252:309);% >5 steps     1.9869    0.0584 (289)
i=8;
[info{i}, files{i}] = organize(7,'male', 'Berlin K', 1.8888, 180317, 1:116);% >5 steps     1.8888    0.0353 (59)
i=9;
[info{i}, files{i}] = organize(8,'male', 'Berlin K', 1.9272, 180317, 117:227);% >5 steps   1.9272    0.0391 (201)
i=10;
[info{i}, files{i}] = organize(8,'male', 'Berlin K', 1.9272, 180318, 1:32); %      1.9190    0.0270 (21)
i=11;
[info{i}, files{i}] = organize(9,'female', 'Berlin K', 1.8881, 180318, 33:58); %Too little data  1.8881    0.0465 (33)
i=12; %No A
[info{i}, files{i}] = organize(10,'female', 'W1118', NaN, 180321, 1:158);% No data
i=13;
[info{i}, files{i}] = organize(11,'female', 'W1118', 2.0420, 180321, 159:286);% >5 steps  2.0420    0.0707 (247)
i=14;
[info{i}, files{i}] = organize(11,'female', 'W1118', 2.0420, 180322, 1:174); %	2.1302    0.0615 (6)
i=15; %No A
[info{i}, files{i}] = organize(12,'female', 'W1118', 2.0354, 180322, 175:489); %needs to be analyzed	2.0354    0.0517 (200) 
i=16;
[info{i}, files{i}] = organize(13,'female', 'Oregon RC', 1.8774, 180324, 1:84); %needs to be confirmed   1.8774    0.0639 (7)
i=17; %No A
[info{i}, files{i}] = organize(14,'female', 'Oregon RC', 2.1137, 180325, 1:112); %Not Analyzed. Bad data. Do not analyze   2.1137    0.1252 (20)
i=18; %No A
[info{i}, files{i}] = organize(15,'female', 'Tac2', NaN, 180614, 1:166); %Data discarded due to tracking error
i=19;
[info{i}, files{i}] = organize(16,'female', 'Tac2', 1.9546, 180614, 167:241); % 1.9546    0.0823 (182)
i=20; %No A
[info{i}, files{i}] = organize(17,'female', 'Tac2', 2.0738, 180625, 1:74); %Leg seems to be cut off. Do not use this data. 2.0738    0.0747 (18)
i=21;
[info{i}, files{i}] = organize(18,'male', 'Berlin K', 2.0646, 180923, 1:232); % 2.0646    0.0912 (32)
i=22;
[info{i}, files{i}] = organize(19,'male', 'Berlin K', 2.0404, 180924, 1:29); % 2.0404    0.0828 (19)
%% find a matching fly and return the info.
for j=1:i
    if ismember([date, vidNum],files{j}, 'rows')
        gender=info{j}{1};
        strain=info{j}{2};
        if strcmp(strain,'W1118')
            if strcmp(gender,'male')
                weight=0.7500;
            elseif strcmp(gender,'female')
                weight=1.1230;
            end 
        elseif strcmp(strain,'Berlin K')
            if strcmp(gender,'male')
                weight=0.6798;
            elseif strcmp(gender,'female')
                weight=1.0809;
            end 
        elseif strcmp(strain,'mutant')
            if strcmp(gender,'male')
                weight=0.6224;
            elseif strcmp(gender,'female')
                weight=0.9687;
            end             
        else
            if strcmp(gender,'male')
                weight=0.7;
            elseif strcmp(gender,'female')
                weight=1.1;
            end             
        end
        
        legLength=info{j}{3};
        %
        flynum=info{j}{4};
        %
        break
    end 
end
end

function [info, files] = organize(flynum, gender, strain, legLength, expDate, dataRange)
info = {gender, strain, legLength, flynum};
files = [ones(length(dataRange),1)*expDate, dataRange'];
end