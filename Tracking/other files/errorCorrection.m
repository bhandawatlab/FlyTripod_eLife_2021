

pathToMatlab = strcat(extractBefore(pwd,'MATLAB'),'MATLAB');
addpath(genpath(pathToMatlab));

pathToAnalysedData = strcat(extractBefore(pwd,'MATLAB'),'Analysis');
d=pathToAnalysedData;
%Creating '\**\shot.mat' for windows and '/**/shot.mat' for matlab
findThis = [filesep '**' filesep 'shot.mat'];
%shotDir = dir(strcat(d,'\**\shot.mat'));
shotDir = dir(strcat(d,findThis));

for i = 1:length(shotDir)
    
    shotName = strcat(shotDir(i).folder,filesep,shotDir(i).name);
    load(shotName);

    if isfield(shot.com.side,'oldError')
        continue
    end
    
    if isfield(shot.com.side,'markers')

        vidFolder = shotDir(i).folder;
        shotPath = strcat(shotDir(i).folder,'\',shotDir(i).name);
        vidDir = strcat(erase(shotPath, '\shot.mat'),'.mat');
        vidDir = strrep(vidDir, 'Analysis', 'Data');
        vidDir = strrep(vidDir,'C:\Users\cc583\Dropbox (Duke Bio_Ea)\Bhandawat Lab\Chun', 'E:\CHUN');
        M = load(vidDir);
        M.data.frames = M.data.frames(:,:,:,shot.frame);
        M.data.timeStamp = M.data.timeStamp(shot.frame);
        
        %vidStartFrame = shot.frame(1);
        sideCoM = [shot.com.bottom.location(:,1), shot.com.worldinfo.floorY - shot.com.side.location(:,2)];
        bottomCoM = shot.com.bottom.location;
        
        errorSide = zeros(size(sideCoM,1),1);
        errorBottom = zeros(size(bottomCoM,1),1);
        for f = 2:size(sideCoM,1)-1 %I am subtracting 1 because shot.com.side/bottom.markers did not contain the data for the markers in the last frame.
            comFuture = sideCoM(f,:);
            markerFuture = shot.com.side.markers{f}';
            comTrackBack=backTraceCoM(M,f,comFuture,markerFuture);
            errorSide(f)=sqrt(sum((sideCoM(f-1,:)-comTrackBack).^2));
            
            if isnan(comTrackBack)
                shotName
                f
            end
            
            comFuture = bottomCoM(f,:);
            markerFuture = shot.com.bottom.markers{f}';
            comTrackBack=backTraceCoM(M,f,comFuture,markerFuture);
            errorBottom(f)=sqrt(sum((bottomCoM(f-1,:)-comTrackBack).^2));
            
            if isnan(comTrackBack)
                shotName
                f
            end
            
        end
        errorSide(1) = errorSide(2);
        errorBottom(1) = errorBottom(2);
        errorSide(end) = errorSide(end-1);
        errorBottom(end) = errorBottom(end-1);
        
        shot.com.side.oldError = shot.com.side.error;
        shot.com.side.oldUncty = shot.com.side.uncty;
        shot.com.side.error = errorSide;
        shot.com.side = rmfield(shot.com.side,'uncty');
        
        shot.com.bottom.oldError = shot.com.bottom.error;
        shot.com.bottom.oldUncty = shot.com.bottom.uncty;
        shot.com.bottom.error = errorBottom;
        shot.com.bottom = rmfield(shot.com.bottom,'uncty');
        
         save(shotName,'shot')
        
    end
    disp([num2str(i/length(shotDir)*100) '%'])
end


function comTrackBack=backTraceCoM(M,fframe,comFuture,markerFuture)

objectFrame = M.data.frames(:,:,:,fframe);%now at frBegin(th) frame

tracker = vision.PointTracker('NumPyramidLevels',3,'MaxBidirectionalError',1);

%COM side view
points = markerFuture;
initialize(tracker,points,objectFrame);

oldPoints = points;

cframe = fframe-1;

frame = M.data.frames(:,:,:,cframe);%-background; %now at frBegin+1

[points, isFound] = step(tracker,frame);

visiblePoints = points(isFound,:);
oldInliers = oldPoints(isFound,:);

if size(visiblePoints, 1) >= 2 % need at least 2 points
    % Estimate the geomatric transformation between the old points
    ...and the new points and eliminate outliers
        [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
        oldInliers, visiblePoints, 'similarity', 'MaxDistance', 1);
    % Apply the transformation to the COM!
    comTrackBack = transformPointsForward(...
        xform, comFuture);
else
    comTrackBack = NaN;
    disp('Incomplete');
end

clear comTrack errr oriPoints markers
release(tracker);

end