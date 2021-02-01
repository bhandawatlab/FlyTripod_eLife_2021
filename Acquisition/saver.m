%This function is used for saving the acquired video into drive.
% @Chanwoo Chun, <cc2465@cornell>

function saver(data)

matdir = pwd;

viddir = ['..' filesep '..' filesep 'Data' filesep 'ProcessedData' filesep datestr(date, 'yymmdd'),'Sideview'];
mkdir(viddir);
cd(viddir);

%///File Naming///&

%Getting rid of the dots that are included with dir() array...
d = dir();
d = d(~ismember({d.name},{'.','..'}));
[m,~] = size(d);
filename=zeros(1,m);
for in = 1:m
    files = d(in).name;
    [~, fName, ~] = fileparts(files);
    %The following will only work when files are simply named with integers
    filename(in)=str2num(fName);
end

%Give it the best name in the world
if isempty(filename)
       newname = 1;
else
       newname=max(filename)+1;
end

disp('Worthy! Saving .mat...');
newnamemat = strcat(num2str(newname), '.mat');
save(newnamemat,'data');
disp('Saved! Resume monitoring... Press any key to stop.');
clear frames;

cd(matdir);
