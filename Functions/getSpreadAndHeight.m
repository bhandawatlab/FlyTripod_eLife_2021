% This function finds L and midstance height of a fly in a given step.
% "data" has an information on the given step, and "shot" has an
% informative on legs.
%
% @Chanwoo Chun, <cc2465@cornell.edu>

function [L, yMid] = getSpreadAndHeight(data,shot)

if ~isfield(shot.leg,'rawPositions')
    L = NaN;
    yMid = NaN;
    return
end
%For entire stance fits
pureTriStarts=data.pureTriStarts;%data.pureTriStarts;%1;
pureTriEnds=data.pureTriEnds-1;%data.pureTriEnds-1;%length(data.time);
rawLegs = [shot.leg.rawPositions.R1';shot.leg.rawPositions.L2.';shot.leg.rawPositions.R3.';...
    shot.leg.rawPositions.L1.';shot.leg.rawPositions.R2.';shot.leg.rawPositions.L3.'];
rawLegs = rawLegs(:,data.source.startFrame:data.source.endFrame);
rawLegs = rawLegs(:,pureTriStarts:pureTriEnds);

eraser = prod(rawLegs,2);
rawLegs = rawLegs(eraser~=0,:);
rawLegs = mean(rawLegs,2);
rawLegsPos = reshape(rawLegs,2,[])';
if size(rawLegsPos,1)~=3
    L = NaN;
    yMid = NaN;
    return
end

%Get anchoring point
comLocation = shot.com.bottom.location;%(data{i}.source.startFrame:data{i}.source.endFrame,:);
comLocation(:,2)=comLocation(:,2)-shot.com.worldinfo.floorY;

bottomLoc = comLocation(data.source.startFrame:data.source.endFrame,:);

allPos=rawLegsPos;

bottomLoc(:,1) = bottomLoc(:,1) - allPos(1,1);
bottomLoc(:,2) = bottomLoc(:,2) - allPos(1,2);
allPos(:,1) = allPos(:,1) - allPos(1,1);
allPos(:,2) = allPos(:,2) - allPos(1,2);

overallAngle = atan2(allPos(3,2),allPos(3,1));

theta = -overallAngle;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
allPos = (R*allPos')';
bottomLoc = (R*bottomLoc')';

%invert correction
if allPos(2,2)<0
    allPos(2,2)=-1*allPos(2,2);
    bottomLoc(:,2)=-1*bottomLoc(:,2);
end
bottomLoc = bottomLoc*24/1984;
allPos=allPos*24/1984;

L = allPos(3,1)/2;

yseries = data.com(:,2);

[~,index] = min(abs(bottomLoc(:,1)-L));
comy = yseries(index);


yMid = comy;


end
