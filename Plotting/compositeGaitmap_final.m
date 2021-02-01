%This code generates composite gait map sorted by speed.
%Segmentation by cycle (NOT segmentation by tripod!)
%
%This code does the following:
% 1) Get cycle start frame and end frame number.
% 2) Plot Figure 3A
%
%%@Chanwoo Chun, Jan. 31, 2021, <cc2465@cornell.edu>

addpath(genpath(['..' filesep '..' filesep '..' filesep 'FlyLocomotion']))
shotDir= dir(['..' filesep '..' filesep '**' filesep 'shot.mat']);

%counter for tripod cycle
cycle1=1;
%counter for nontripod cycle
cycle2=-1;
%counter for all cycles
cycleNum=1;
%cycleNum acts as an index for cycleCell
cycleMat=nan(6,2000);
flyNums=nan(1,2000);
compGaitCell=cell(2,10000);

cycleCount = nan(length(shotDir),1);
ddd=0;
for i = 1:length(shotDir)
    cyclePerVid=0;
    shotName = strcat(shotDir(i).folder,filesep,shotDir(i).name);
    load(shotName);
    
    t = shot.timeStamp';
    dxdt=shot.com.bottom.dxdt*24/1984;
    dydt=shot.com.bottom.dydt*24/1984;
    comVel = (dxdt.^2+dydt.^2).^(1/2);
    mmCoMVel = movmean(comVel,10); %moving mean of CoM velocity.
    
    [gender, strain, weight, legLength, flynum] = getFlyInfo(shotDir(i).folder);
    
    valid=shot.validity.valid;
    legs = [shot.leg.R1.';shot.leg.R2.';shot.leg.R3.';shot.leg.L1.';shot.leg.L2.';shot.leg.L3.'];
    
    %The below four lines were needed because of ONE strange data we had.
    %In the data, one leg was tracked outside the valid portion of the
    %data. This messed up the analysis. Therefore, the data is omitted.
    summed=sum(legs,1);
    if any(summed(logical(1-valid)))
        ddd=ddd+1
        continue
    end
    
    [vStarts, vEnds] = getStartsAndEnds(valid);
    
    SE = cell(12,length(vStarts));
    for legNum = 1:6
        [starts, ends] = getStartsAndEnds(legs(legNum,:));
        [starts, ends] = removeInvalidStartsAndEnds(starts,ends,vStarts,vEnds,'option2');
        
        for vi=1:length(vStarts)
            %SE's 1st~6th cells: starting indices
            %SE's 7th~12th cells: ending indices
            SE{legNum,vi}=starts(starts>=vStarts(vi)&starts<=vEnds(vi));
            SE{legNum+6,vi}=ends(ends>=vStarts(vi)&ends<=vEnds(vi));
        end
    end
    
    %Loop through valid sections
    for vi=1:length(vStarts)
        
        %SOI: starts of interest.
        SOI=SE(1:6,vi);
        EOI=SE(7:12,vi);
        infoCell={1,6};
        for j=1:6
            if any([1 3 5]==j)
                leg_triMark=1;
            else
                leg_triMark=-1;
            end
            
            %loi means leg of interest
            soi_loi=SOI{j};
            eoi_loi=EOI{j};
            infoCell{j}=[j*ones(1,length(soi_loi)); soi_loi; eoi_loi; leg_triMark*ones(1,length(soi_loi))];
            
            if isempty(infoCell{j})
                infoCell{j}=[[];[];[]];
            end
        end
        info = cell2mat(infoCell);
        if size(info,2)<3
            %this is the case where there was less than three stance starts.
            continue
        end
        [~,sortI] = sort(info(2,:));
        sequence = info(:,sortI);
        
        %sequence matrix tells which leg (represented with leg number at
        %the first row) is landed at which frame (second row), sorted by
        %time of incident. 3rd row stance ends. The 4th low indicates which
        %tripod the leg came from.
        %Now, identify if there are tripods.
        temp2=movsum(sequence(4,:),3,'Endpoints','discard');
        logicalI1=temp2==3;
        tripod1Starts=find(logicalI1);
        logicalI2=temp2==-3;
        tripod2Starts=find(logicalI2);
        allTriStarts = sort([tripod1Starts tripod2Starts]);
        
        %In markTriStarts, the start of right tripod will be denoted as 1,
        %and that of left tripod as -1.
        markTriStarts=[logicalI1-logicalI2 0 0];
        %markTriStarts going in as 5th row.
        sequence = [sequence;markTriStarts];
        
        cycleMarker = zeros(1,size(sequence,2));
        %Loop through all tripod start indices
        for k = 1:length(allTriStarts)
            refPoint = allTriStarts(k);
            if cycleMarker(refPoint)~=0
                %In this case, cycle was already identified in this given
                %frame.
                continue
            end
            
            
            %1) identify first leg
            refLeg=sequence(1,refPoint);
            %2) get frame number when the same leg comes again.
            sameLegFrames=sequence(2,sequence(1,:)==refLeg);
            endFrame=sameLegFrames(find(sameLegFrames>sequence(2,refPoint),1));
            
            refFrame = sequence(2,refPoint);
            
            stanceDonePoint=[];
            for q = refPoint:size(sequence,2)
                legs_till_now=unique(sequence(1,refPoint:q));
                %If the following is true, all legs have gone through
                %stance start phase.
                if length(legs_till_now)==6
                    stanceDonePoint = q;
                    stanceDoneFrame = sequence(2,stanceDonePoint);
                    break
                end
            end
            
            if isempty(stanceDonePoint)
                continue
            end
            
            %So, endFrame is a point in time when one cycle is completed,
            %and stanceDoneFrame is a point in time when all leg have gone
            %through their stance phases. refFrame is when gait cycle
            %started.
            cycleMarker(refPoint:stanceDonePoint)=cycle1;
            
            compGait=getGaitMap(sequence,refPoint,stanceDonePoint,t);
            
            %get average speed over this time period.
            avgSpeed=mean(comVel(refFrame:stanceDoneFrame));
            
            compGaitCell{1,cycleNum}=avgSpeed;
            compGaitCell{2,cycleNum}=compGait;
            
            flyNums(cycleNum) = flynum;
            
            cycle1=cycle1+1;
            cycleNum=cycleNum+1;
            if ~isempty(endFrame)
                cyclePerVid=cyclePerVid+1;
            end
        end
        
        %Now check if there is any instance where no cycle was assigned.
        %If this instance includes all six legs and if the cycle duration
        %can be calculated, we call it unidentified gait, and than assign
        %gait cycle number (cycle2), which is a negative value.
        [starts, ends] = getStartsAndEnds(cycleMarker==0);
        
        %Starts and ends of a time frame where one or more unidentified
        %gaits could possibly be found.
        %i4s: index for starts
        for i4s = 1:length(starts)
            %Scan through all points. Any point can be a gait starting
            %point.
            for u_refPoint = starts(i4s):ends(i4s)
                if cycleMarker(u_refPoint)~=0
                    %In this case, a cycle was already identified in this 
                    %given frame.
                    continue
                end
                %Now check if there is a complete gait and then identify
                %it (mark it).
                u_stanceDonePoint=[];
                for q = u_refPoint:ends(i4s)
                    legs_till_now=unique(sequence(1,u_refPoint:q));
                    %If the following is true, all legs have gone through
                    %stance start phase.
                    if length(legs_till_now)==6
                        u_stanceDonePoint = q;
                        u_stanceDoneFrame = sequence(2,u_stanceDonePoint);
                        break
                    end
                end
                if isempty(u_stanceDonePoint)
                    continue
                end
                
                
                %1) identify first leg
                u_refLeg=sequence(1,u_refPoint);
                
                u_sameLegFrames=sequence(2,sequence(1,:)==u_refLeg);
                u_endFrame=u_sameLegFrames(find(u_sameLegFrames>sequence(2,u_refPoint),1));
                
                u_refFrame = sequence(2,u_refPoint);
                
                %Now, if this iteration made it to this point, offically
                %mark it as a complete unidentified gait.
                cycleMarker(u_refPoint:u_stanceDonePoint)=cycle2;
                
                u_compGait=getGaitMap(sequence,u_refPoint,u_stanceDonePoint,t);
                
                %get average speed over this time period.
                u_avgSpeed=mean(comVel(u_refFrame:u_stanceDoneFrame));
                
                compGaitCell{1,cycleNum}=u_avgSpeed;
                compGaitCell{2,cycleNum}=u_compGait;
                
                flyNums(cycleNum) = flynum;
                
                
                cycle2=cycle2-1;
                cycleNum=cycleNum+1;
                if ~isempty(u_endFrame)
                    cyclePerVid=cyclePerVid+1;
                end
            end
        end
        
        %cycleMarker going in as 5th row
        sequence = [sequence;cycleMarker];
        
    end
    
    cycleCount(i) = cyclePerVid;
end

compGaitCell(:,cycleNum:end)=[];
[~,compGaitSortIdx]=sort(cell2mat(compGaitCell(1,:)));
compGaitCell=compGaitCell(:,compGaitSortIdx);

slide = 1/size(compGaitCell,2)*0.9;
figure
hold on
for g = 1:size(compGaitCell,2)
    gaitMap = compGaitCell{2,g}(1:6,:);
    timeStamp = compGaitCell{2,g}(7,:);
    gaitMap(1,:) = gaitMap(1,:)*6;
    gaitMap(5,:) = gaitMap(5,:)*5;
    gaitMap(3,:) = gaitMap(3,:)*4;
    gaitMap(4,:) = gaitMap(4,:)*3;
    gaitMap(2,:) = gaitMap(2,:)*2;
    gaitMap(6,:) = gaitMap(6,:)*1;
    
    gaitMap = gaitMap - g*slide;
    handl = plot(timeStamp,gaitMap);
    set(handl, {'color'}, {[242 110 37]/255; [65 190 238]/255; [247 212 189]/255;[29 147 192]/255; [247 168 122]/255; [148 215 242]/255});
end
hold off
xlim([-0.2 0.4])
set(gcf,'renderer','painters');

flyNumCount=flyCount(flyNums);

numVidsAtLeastThreeSteps = numel(find(cycleCount>=3));

disp(['#vids with >=3 steps: ' num2str(numVidsAtLeastThreeSteps)])
disp(['#vids: ' num2str(length(cycleCount))])

function compGait=getGaitMap(sequence,refPoint,stanceDonePoint,t)
%Make composite gait map data
%lastFrame is a frame number when the last leg's stance ends.
refLeg = sequence(1,refPoint);

refFrame=sequence(2,refPoint);
lastFrame=max(sequence(3,refPoint:stanceDonePoint));
numFrame=lastFrame-refFrame+1;
compGait = NaN(7,numFrame);

timeSeries=t(refFrame:lastFrame);
compGait(7,:)=timeSeries;
for p = refPoint:stanceDonePoint
    %ln: leg number
    ln = sequence(1,p);
    sStart = sequence(2,p)-refFrame+1;
    sEnd = sequence(3,p)-refFrame+1;
    compGait(ln,sStart:sEnd)=1;
    

end

if sequence(4,refPoint) == -1
    compGait = compGait([4 5 6 1 2 3 7],:);
end

proLeg = compGait(1,:);
proLeg(isnan(proLeg))=0;
[proStart, ~] = getStartsAndEnds(proLeg);
proStart=min(proStart);
compGait(7,:)=compGait(7,:)-compGait(7,proStart);

end