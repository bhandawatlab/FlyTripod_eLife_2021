%This is the code that will fit ARSLIP, SLIP and IP models.
%Before fitting, the all data are segmented into tripod stance phases.
%The for-loop starting with "for i = 1:length(shotDir)" takes care of the
%segmentation.
%
%@Chanwoo Chun, <cc2465@cornell.edu>

clear all
close all

addpath(genpath(['..' filesep '..' filesep '..' filesep 'FlyLocomotion']))
shotDir= dir(['..' filesep '..' filesep '**' filesep 'shot.mat']);

allVid = cell(1,1000);
sequenceCell = cell(1,100000);
idx = 1;
aIdx = 1;
for i = 1:length(shotDir)
%     if isBadData(shotDir(i).folder,'velPlot')
%         continue
%     end
    
    shotName = strcat(shotDir(i).folder,filesep,shotDir(i).name);
    load(shotName);
    
    valid=shot.validity.valid;
    validLegs = [shot.leg.R1.';shot.leg.R2.';shot.leg.R3.';shot.leg.L1.';shot.leg.L2.';shot.leg.L3.'];
    time = shot.timeStamp';%*1000;
    dx = gradient(shot.com.bottom.location(:,1));
    dy = gradient(shot.com.bottom.location(:,2));
    dz = gradient(shot.com.side.location(:,2));
    dt = gradient(shot.timeStamp);
    dxdt = dx./dt*24/1984;
    dydt = dy./dt*24/1984;
    dzdt = dz./dt*24/1984;
    comVel = (dxdt.^2+dydt.^2).^(1/2);
    comVel3D = [comVel dzdt];%(dxdt.^2+dydt.^2+dzdt.^2).^(1/2);
    mmCoMVel = movmean(comVel,10); %moving mean of CoM velocity.
    
    %calculate distance traveled
    dx=[0; diff(shot.com.bottom.location(:,1))];
    dy=[0; diff(shot.com.bottom.location(:,2))];
    dtravel = (dx.^2+dy.^2).^(1/2);
    dista = cumsum(dtravel)';
    distamm = dista.*24/1984;
    
    boundaryTemp = cell(1,2);
    for flip = 1:2
        if flip == 2
            validLegs = fliplr(validLegs);
            time = fliplr(time);
            valid = fliplr(valid);
            comVel = fliplr(comVel);
            mmCoMVel = fliplr(mmCoMVel);
        end
        %determin starts and ends of the valid portion
        v=valid';
        w = [false v~=0 false]; %// "close" v with zeros, and transform to logical
        validStarts = find(w(2:end) & ~w(1:end-1)); %// find starts of runs of non-zeros
        validEnds = find(~w(2:end) & w(1:end-1))-1; %// find ends of runs of non-zeros
        
        %if length(validStarts)<2
        %    continue
        %end
        
        %sequence matrix
        stanceMarked = zeros(6,length(valid),2);
        swingMarked = zeros(6,length(valid));
        for h = [1 5 3 4 2 6]
            
            v = validLegs(h,:); %// data <<<<
            w = [false v~=0 false]; %// "close" v with zeros, and transform to logical
            starts = find(w(2:end) & ~w(1:end-1)); %// find starts of runs of non-zeros
            ends = find(~w(2:end) & w(1:end-1))-1; %// find ends of runs of non-zeros
            
            startsDur = starts;
            endsDur = ends;
            
            %The following statements crop out the first and last step if they
            ...are incomplete.
                [startsDur, iS] = setdiff(startsDur,validStarts); %remove starting indices that are the same as valid start indices
            endsDur = endsDur(iS); %remove corresponding ending indices FOR DURATION CALCULATION
            [endsDur,iE] = setdiff(endsDur,validEnds); %same for the ending indices
            startsDur = startsDur(iE); %same FOR DURATION CALCULATION
            
            stanceDuration = time(endsDur)-time(startsDur);
            
            %Now that we know the stance duration, find the cycle duration and
            %get duty factor.
            dutyFactor = zeros(1,length(startsDur));
            for sdi = 1:length(startsDur)
                strideEndIdx =  min(starts(starts>startsDur(sdi)))-1;
                
                if isempty(strideEndIdx)
                    strideStartIdx = max(ends(ends<endsDur(sdi)))+1;
                    
                    if isempty(strideStartIdx)
                        continue
                    end
                    
                    dutyFactor(sdi) = stanceDuration(sdi)/(time(endsDur(sdi)) - time(strideStartIdx));
                    continue
                end
                
                dutyFactor(sdi) = stanceDuration(sdi)/(time(strideEndIdx)-time(startsDur(sdi)));
            end
            ends = ends + 1; %Make theses indices mean the swing start frame.
            
            switch h
                case 1
                    label = 1;
                case 2
                    label = -2;
                case 3
                    label = 3;
                case 4
                    label = -1;
                case 5
                    label = 2;
                case 6
                    label = -3;
            end
            
            stanceMarked(h,starts,1) = label;
            stanceMarked(h,startsDur,2) = stanceDuration;%dutyFactor;
            swingMarked(h,ends) = label;
        end
        
        
        %look into each valid section
        boundariesAll = cell(1,length(validStarts));
        for vi = 1:length(validStarts)
            %Swing does not seem to work correctly.
            rightTripods = cell(1,1000);
            rtIdx = 1;
            leftTripods = cell(1,1000);
            ltIdx = 1;
            markedCropped =  stanceMarked(:,validStarts(vi):validEnds(vi),:);%%%%%%%%%%%%%%%%%STANCE OR SWING?%%%%%%%%%%%%
            timeCropped = time(validStarts(vi):validEnds(vi));
            frameCell = cell(1,size(markedCropped,2));
            for f = 1:size(markedCropped,2)
                framNum = validStarts(vi)+f-1;
                frameMark=squeeze(markedCropped(:,f,:));
                frameMark(frameMark(:,1)==0,:)=[];
                frameDF = frameMark(:,2)';
                frameMark = frameMark(:,1)';
                
                if isequal(frameMark,[-1 3])
                    frameMark = [3 -1]; %frameMark = fliplr(frameMark);
                    frameDF = fliplr(frameDF);
                end
                if isequal(frameMark,[1 -3])
                    frameMark = [-3 1];
                    frameDF = fliplr(frameDF);
                end
                
                frameTime = timeCropped(f);
                frameTimes = repmat(frameTime,1,length(frameMark));
                frameVel = mmCoMVel(f);
                frameVels = repmat(frameVel,1,length(frameMark));
                framNum = repmat(framNum,1,length(frameMark));
                
                frameCell{f} = [frameMark; frameTimes; framNum; frameDF];
            end
            sequenceMat = cell2mat(frameCell);
            
            if isempty(sequenceMat)
                continue
            end
            
            tripods = {[1 2 3],[1 3 2],[2 1 3], [2 3 1], [3 1 2], [3 2 1]};
            triMark = zeros(1,size(sequenceMat,2));
            for j = 1:6
                for rightleft = [1 -1]
                    trIdx = strfind(sequenceMat(1,:),tripods{j}*rightleft);
                    
                    for k = 1:length(trIdx)
                        intIdx = [trIdx(k),trIdx(k)+1,trIdx(k)+2]; %Interested Idex (index of interest)
                        if any(triMark(intIdx))
                            continue
                        end
                        triMark(intIdx)=1;
                        
                        frames = sequenceMat(3,intIdx);
                        legNum=sequenceMat(1,intIdx);
                        legNumTemp=legNum;
                        legNum(legNumTemp==-2)=2;
                        legNum(legNumTemp==-1)=4;
                        legNum(legNumTemp==2)=5;
                        legNum(legNumTemp==-3)=6;
                        
                        triGroup = zeros(3,size(validLegs,2));
                        good = true;
                        for p = 1:3
                            validLegLabel = bwlabel(validLegs(legNum(p),:));
                            thisStance = validLegLabel==validLegLabel(frames(p));
                            triGroup(p,:) = thisStance;
                        end
                        
                        if rightleft == 1
                            rightTripods{rtIdx} = sum(triGroup,1);
                            rtIdx = rtIdx + 1;
                        else
                            leftTripods{ltIdx} = sum(triGroup,1);
                            ltIdx = ltIdx + 1;
                        end
                    end
                end
            end
            
            rightTripods(rtIdx:end) = [];
            leftTripods(ltIdx:end) = [];
            boundaries = zeros(1,1000);
            brIdx = 1;
            for ri = 1:rtIdx-1
                RT = rightTripods{ri}>=1;
                for li = 1:ltIdx-1
                    LT = leftTripods{li}>=1;
                    intersection = RT.*LT;
                    %Remove the first and last incomplete tripod boundaries
                    if ~isempty(intersect(find(intersection),[validStarts validEnds]))
                        continue
                    end
                    midpoint = round((find(intersection,1,'last')+find(intersection,1))/2);
                    if ~isempty(midpoint)
                        boundaries(brIdx) = midpoint;
                        brIdx = brIdx+1;
                    end
                end
            end
            boundaries(brIdx:end)=[];
            boundariesAll{vi} = boundaries;
        end
        
        boundariesAllMat=cell2mat(boundariesAll);
        if isempty(boundariesAllMat)
            boundariesAllMat=double.empty(1,0);
        end
        [~,order]=sort(boundariesAllMat(1,:));
        boundariesAllMat = boundariesAllMat(:,order);
        
        if flip == 2
            boundariesAllMat = size(validLegs,2)-boundariesAllMat+2;
            validLegs = fliplr(validLegs);
            time = fliplr(time);
            valid = fliplr(valid);
            comVel = fliplr(comVel);
            mmCoMVel = fliplr(mmCoMVel);
        end
        boundaryTemp{flip} = boundariesAllMat;
    end
    
    bori=boundaryTemp{1};
    bflip=boundaryTemp{2};
    
    if ~isempty(bori)
        bflip = bflip(:,bflip(1,:)<(min(bori(1,:))-3));
    else
        disp('whelp')
    end
    boundary = [bflip bori];
    [~,order]=sort(boundary(1,:));
    boundary = boundary(:,order);
    
    %% Now pack up for fitting.
    markedLegs = [shot.leg.R1.';shot.leg.L2.';shot.leg.R3.';shot.leg.L1.'*10;shot.leg.R2.'*10;shot.leg.L3.'*10];
    stanceBounds=boundary(1,:);
    OVidx=1;
    oneVidCell = cell(1,length(stanceBounds));
    
    if length(stanceBounds)==1
        continue
    end
    
    for fi=1:length(stanceBounds)-1
        if 0 == prod(shot.validity.valid(stanceBounds(fi):stanceBounds(fi+1)))
            continue
        end
        
        s.height = shot.com.side.location(stanceBounds(fi):stanceBounds(fi+1),2)*24/1984;
        s.time = shot.timeStamp(stanceBounds(fi):stanceBounds(fi+1));
        if rem(length(s.time),2)==0
            meanMid = (s.time(length(s.time)/2)+s.time(length(s.time)/2+1))/2;
            s.time = s.time-meanMid;
        else
            s.time = s.time - s.time(round(length(s.time)/2));
        end
        
        s.source.gender=nan;
        s.source.strain=nan;
        s.source.weight=nan;
        s.source.legLength=nan;
        s.source.flynum=nan;
        try
        [s.source.gender, s.source.strain, s.source.weight, legLength, flynum] = getFlyInfo(shotDir(i).folder);
        catch
            disp('source info are nan')
        end
        s.source.legLength = legLength;
        s.source.flynum = flynum;
        
        s.com(:,1) = abs(distamm(stanceBounds(fi):stanceBounds(fi+1))-distamm(stanceBounds(fi)));
        s.com(:,2) = s.height; %shot.com.side.location(starts(q):ends(q),2).*24/1984;
        
        [filepathTemp,vidNum,~] = fileparts(shotDir(i).folder);
        [~,folder,~]=fileparts(filepathTemp);
        s.source.videoDir = [folder filesep vidNum];
        
        s.source.startFrame = stanceBounds(fi);
        s.source.endFrame = stanceBounds(fi+1);
        avgVelStance=mean(comVel((stanceBounds(fi):stanceBounds(fi+1))));
        s.avgVel = avgVelStance;
        s.sideErrormm = shot.com.side.error(stanceBounds(fi):stanceBounds(fi+1))*24/1984;
        s.bottomErrormm = shot.com.bottom.error(stanceBounds(fi):stanceBounds(fi+1))*24/1984;
        s.source.bodyLength = shot.com.worldinfo.bodyLength;
        
        s.vel = comVel(stanceBounds(fi):stanceBounds(fi+1)); %mm/s
        s.vel3D = comVel3D(stanceBounds(fi):stanceBounds(fi+1),:);
 
        normalTime = (s.time - s.time(1))*2/(s.time(end)-s.time(1))-1;
        
        s.theTime = -1:0.01:1;

        s.height_interp = interp1(normalTime,s.height,s.theTime);
        s.vel_interp = interp1(normalTime,s.vel,s.theTime);
        
        markedLegsTemp = markedLegs(:,stanceBounds(fi):stanceBounds(fi+1));
        triIndicator = sum(markedLegsTemp,1);
        v= triIndicator==3 | triIndicator==30;
        w = [false v~=0 false]; %// "close" v with zeros, and transform to logical
        pureTriStarts = find(w(2:end) & ~w(1:end-1)); %// find starts of runs of non-zeros
        pureTriEnds = find(~w(2:end) & w(1:end-1)); %// find ends of runs of non-zeros
        s.pureTriStarts = pureTriStarts;
        s.pureTriEnds = pureTriEnds;
        
        s.PTR = abs(time(pureTriStarts)-time(pureTriEnds))/abs(s.time(end)-s.time(1))*100;
        
        s.height_interp = interp1(normalTime,s.height,s.theTime);
        s.vel_interp = interp1(normalTime,s.vel,s.theTime);
        
        if isfield(shot.leg,'rawPositions')
        s.leg.rawLegPos.R1 = shot.leg.rawPositions.R1(stanceBounds(fi):stanceBounds(fi+1),:);
        s.leg.rawLegPos.R2 = shot.leg.rawPositions.R2(stanceBounds(fi):stanceBounds(fi+1),:);
        s.leg.rawLegPos.R3 = shot.leg.rawPositions.R3(stanceBounds(fi):stanceBounds(fi+1),:);
        s.leg.rawLegPos.L1 = shot.leg.rawPositions.L1(stanceBounds(fi):stanceBounds(fi+1),:);
        s.leg.rawLegPos.L2 = shot.leg.rawPositions.L2(stanceBounds(fi):stanceBounds(fi+1),:);
        s.leg.rawLegPos.L3 = shot.leg.rawPositions.L3(stanceBounds(fi):stanceBounds(fi+1),:);
        s.leg.R1 = shot.leg.R1(stanceBounds(fi):stanceBounds(fi+1),:);
        s.leg.R2 = shot.leg.R2(stanceBounds(fi):stanceBounds(fi+1),:);
        s.leg.R3 = shot.leg.R3(stanceBounds(fi):stanceBounds(fi+1),:);
        s.leg.L1 = shot.leg.L1(stanceBounds(fi):stanceBounds(fi+1),:);
        s.leg.L2 = shot.leg.L2(stanceBounds(fi):stanceBounds(fi+1),:);
        s.leg.L3 = shot.leg.L3(stanceBounds(fi):stanceBounds(fi+1),:); 
        end
        
        [s.tripod.L, s.tripod.yMid] = getSpreadAndHeight(s,shot);
        
        
        if numel(s.PTR) ~= 1
            s.PTR = 0;
        end
        
        oneVidCell{OVidx}=s;
        
        OVidx=OVidx+1;
        clear s
    end
    oneVidCell(:,OVidx:end)=[];
    allVid{aIdx}=oneVidCell;
    
    aIdx=aIdx+1;
    
    stanceBoundsTimes=time(boundary(1,:));
end

allVid(:,aIdx:end)=[];
sequenceCell(:,idx:end) = [];

stepData=[allVid{:}];

%% Fit models and save.
filepath = ['..' filesep '..' filesep 'Data' filesep 'StepData' filesep 'stepData1500.mat'];

parpool('local',4)

stepData = fitLooper(stepData,'ARSLIP');
save(filepath,'stepData')

stepData = fitLooper(stepData,'SLIP');
save(filepath,'stepData')

stepData = fitLooper(stepData,'IP');
save(filepath,'stepData')

delete(gcp);


function new_stepData = fitLooper(stepData,modelName)
PTRthres=25;
cnt=0;
for i = 1:length(stepData)
    if isempty(stepData{i}.pureTriStarts)
        continue
    end
    if stepData{i}.PTR<PTRthres || isempty(stepData{i}.PTR)
        continue
    end

    cnt = cnt+1;

end
cnt
ppm = ParforProgMon('progress', cnt);
parfor st = 1:length(stepData)
    stepFit = stepData{1,st};
    if isempty(stepData{1,st}.pureTriStarts)
        continue
    end
    if stepFit.PTR<PTRthres || isempty(stepFit.PTR)
        continue
    end
    
    switch modelName
        case 'ARSLIP'
            stepFit = fitARSLIP(stepFit,st);
        case 'SLIP'
            stepFit = fitSLIP(stepFit,st);
        case 'IP'
            stepFit = fitIP(stepFit,st);
    end
    
    stepData{1,st} = stepFit;
    
    ppm.increment();
end
new_stepData=stepData;
end


