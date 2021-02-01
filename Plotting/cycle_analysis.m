%This code performs gait analysis: time and phase based delay between legs
%Fig 3B,C,E,F
%
%Segmentation by cycle (NOT segmentation by tripod!)
%
%@Chanwoo Chun, Jan. 31, 2021, <cc2465@cornell.edu>

global cycleNum

addpath(genpath(['..' filesep '..' filesep '..' filesep 'FlyLocomotion']))
shotDir= dir(['..' filesep '..' filesep '**' filesep 'shot.mat']);


% For each cycle, the legs' stance start frame numbers will be recorded
% into a single matrix (e.g. A).
% A = [0.3  1.1 0.4 1.2 0.6 1.4] <-stance start times
%       R1   R2  R3  L1  L2  L3  <-of these legs.

bound = [0 sqrt(sum((2*[pi pi pi pi]).^2))];
%
relativeExp=cell(1,10000);
rNum = 1;
relaPAvg = zeros(5,10000);
%counter for tripod cycle
cycle1=1;
%counter for nontripod cycle
cycle2=-1;
%counter for all cycles
cycleNum=1; 
%cycleNum acts as an index for cycleCell
cycleMat=nan(6,2000);

deltaMatTime=nan(4,2000);
deltaMatPhase=nan(4,2000);

distTime = nan(5,2000);
distPhase = nan(5,2000); 
GDI_E_list = nan(4,2000);
distTriTet = nan(2,2000);
deltasToR1Mat = nan(5,2000); 
menIdx = nan(2,2000);

flyNums=nan(1,2000);

for i = 1:length(shotDir)
    disp(num2str(i/length(shotDir)*100))
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
    
    [vStarts, vEnds] = getStartsAndEnds(valid);
    
    analyticalSignals=[];
    if isfield(shot.leg,'legPositionPlot')
        legPos=shot.leg.legPositionPlot;
        %Get analytical signals of the legs
        analyticalSignals=getAnalyticalSignal(legPos,vStarts,vEnds);
    else
        continue
    end
    
    SE = cell(12,length(vStarts));
    for legNum = 1:6
        [starts, ends] = getStartsAndEnds(legs(legNum,:));
        [starts, ends] = removeInvalidStartsAndEnds(starts,ends,vStarts,vEnds,'option1');
        
        for vi=1:length(vStarts)
            %SE's 1st~6th cells: starting indices
            %SE's 7th~12th cells: ending indices
            SE{legNum,vi}=starts(starts>=vStarts(vi)&starts<=vEnds(vi));
            SE{legNum+6,vi}=ends(ends>=vStarts(vi)&ends<=vEnds(vi));
        end
    end
    
    %Loop through valid sections
    for vi=1:length(vStarts)
        expPhases=analyticalSignals(:,vStarts(vi):vEnds(vi));
        if any(~isnan(sum(expPhases,1))) %size(expPhases,2)>2
        [pExpT, ~] = unwrapAbout(expPhases,7,[1 3 5]);
        relaP = [pExpT(2,:)-pExpT(1,:);pExpT(3,:)-pExpT(1,:);pExpT(4,:)-pExpT(1,:);pExpT(5,:)-pExpT(1,:);pExpT(6,:)-pExpT(1,:)];
        end
        rNum=rNum+1;
        
        %SOI: starts of interest.
        SOI=SE(1:6,vi);
        
        infoCell={1,6};
        for j=1:6
            if any([1 3 5]==j)
                leg_triMark=1;
            else
                leg_triMark=-1;
            end
            
            %loi means leg of interest
            soi_loi=SOI{j};
            infoCell{j}=[j*ones(1,length(soi_loi)); soi_loi; leg_triMark*ones(1,length(soi_loi))];
            
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
        %time of incident. The third low indicates which tripod the leg
        %came from.
        %Now, identify if there are tripods.
        temp2=movsum(sequence(3,:),3,'Endpoints','discard');
        logicalI1=temp2==3;
        tripod1Starts=find(logicalI1);
        logicalI2=temp2==-3;
        tripod2Starts=find(logicalI2);
        allTriStarts = sort([tripod1Starts tripod2Starts]);
        
        %In markTriStarts, the start of right tripod will be denoted as 1,
        %and that of left tripod as -1.
        markTriStarts=[logicalI1-logicalI2 0 0];
        %markTriStarts going in as 4th row.
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
            
            %Get cycle duration
            %1) identify first leg
            refLeg=sequence(1,refPoint);
            %2) get frame number when the same leg comes again.
            sameLegFrames=sequence(2,sequence(1,:)==refLeg);
            endFrame=sameLegFrames(find(sameLegFrames>sequence(2,refPoint),1));
            if isempty(endFrame)
                continue
            end
            %3) get cycle duration
            refFrame = sequence(2,refPoint);
            cycleDuration = t(endFrame)-t(refFrame); 
            
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
            
            %get average speed over this time period.
            avgSpeed=mean(comVel(refFrame:endFrame));
            
            
            [deltas0, deltas1, deltasTT, deltasToR1, distT, GDI, GDI_E] = getDataForPlots(sequence,refPoint,stanceDonePoint,cycleDuration,t,avgSpeed);
            
            [tri_md, tet_md]=getMenIdx(legs(:,refFrame:endFrame));
            menIdx(:,cycleNum) = [tri_md; tet_md];
            
            relaPhase=nan;
            phaseGDI=nan;
            distP=nan;
            [relaPhase, phaseGDI, distP] = getPhaseGDI(shot,analyticalSignals,t,refFrame,endFrame,refLeg,avgSpeed);
            
            %plotGaitMapSortByGDI(legs,refFrame,endFrame,distP,t,cycleNum);
            
%             if cycleNum==438 || cycleNum==468
%             end
            %get
            %Pack up
            %cycleCell
            %   row 1: speed
            cycleMat(1,cycleNum) = avgSpeed;
            %   row 2~4: deltas for plotting
            cycleMat(2:4,cycleNum) = deltas0;
            %   row 5: GDI
            cycleMat(5,cycleNum) = GDI;
            %   row 6: phase based GDI
            cycleMat(6,cycleNum) = phaseGDI;
            
            deltaMatTime(:,cycleNum) = deltas1';
            distTime(:,cycleNum) = distT;
            distPhase(:,cycleNum) = distP;
            GDI_E_list(:,cycleNum) = GDI_E;
            distTriTet(:,cycleNum) = deltasTT;
            deltasToR1Mat(:,cycleNum) = deltasToR1;
            
            if any(~isnan(sum(expPhases,1))) %size(expPhases,2)>2
            pExpShiftAvgAll = mean(relaP(:,refFrame-vStarts(vi)+1:endFrame-vStarts(vi)+1),2,'omitnan');
            relaPAvg(:,cycleNum) = pExpShiftAvgAll;%/(2*pi);
            end
            
            flyNums(cycleNum) = flynum;
            
            cycle1=cycle1+1;
            cycleNum=cycleNum+1;
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
                
                %Now, check if I can get gait cycle duration.
                %1) identify first leg
                u_refLeg=sequence(1,u_refPoint);
                %2) get frame number when the same leg comes again.
                u_sameLegFrames=sequence(2,sequence(1,:)==u_refLeg);
                u_endFrame=u_sameLegFrames(find(u_sameLegFrames>sequence(2,u_refPoint),1));
                if isempty(u_endFrame)
                    continue
                end
                %3) get cycle duration
                u_refFrame = sequence(2,u_refPoint);
                u_cycleDuration = t(u_endFrame)-t(u_refFrame);
                
                %Now, if this iteration made it to this point, offically
                %mark it as a complete unidentified gait.
                cycleMarker(u_refPoint:u_stanceDonePoint)=cycle2;
                
                %get average speed over this time period.
                u_avgSpeed=mean(comVel(u_refFrame:u_endFrame));
                
            
                %Now calculate necessary parameters.
                [u_deltas0, u_deltas1, u_deltasTT, u_deltasToR1, u_distT, u_GDI, u_GDI_E] = getDataForPlots(sequence,u_refPoint,u_stanceDonePoint,u_cycleDuration,t,u_avgSpeed);
                
                [u_tri_md, u_tet_md]=getMenIdx(legs(:,u_refFrame:u_endFrame));

                menIdx(:,cycleNum) = [u_tri_md; u_tet_md];
                
                u_relaPhase=nan;
                u_phaseGDI=nan;
                u_distP=nan;
                [u_relaPhase, u_phaseGDI, u_distP] = getPhaseGDI(shot,analyticalSignals,t,u_refFrame,u_endFrame,u_refLeg,u_avgSpeed);
                
                %plotGaitMapSortByGDI(legs,u_refFrame,u_endFrame,u_distP,t,cycleNum);
               
                %Pack up
                %cycleCell
                %   row 1: speed
                cycleMat(1,cycleNum) = u_avgSpeed;
                %   row 2~4: deltas for plotting
                cycleMat(2:4,cycleNum) = u_deltas0;
                %   row 5: time based GDI
                cycleMat(5,cycleNum) = u_GDI;
                %   row 6: phase based GDI
                cycleMat(6,cycleNum) = u_phaseGDI;
                
                deltaMatTime(:,cycleNum) = u_deltas1';
                distTime(:,cycleNum) = u_distT;
                distPhase(:,cycleNum) = u_distP;
                GDI_E_list(:,cycleNum) = u_GDI_E;
                distTriTet(:,cycleNum) = u_deltasTT;
                deltasToR1Mat(:,cycleNum) = u_deltasToR1;
                
                flyNums(cycleNum) = flynum;
                
                if any(~isnan(sum(expPhases,1))) %size(expPhases,2)>2
                u_pExpShiftAvgAll = mean(relaP(:,u_refFrame-vStarts(vi)+1:u_endFrame-vStarts(vi)+1),2,'omitnan');
                relaPAvg(:,cycleNum) = u_pExpShiftAvgAll;%/(2*pi);
                end
                
                cycle2=cycle2-1;
                cycleNum=cycleNum+1;
            end
        end
        
        %cycleMarker going in as 5th row
        sequence = [sequence;cycleMarker];
        
    end
    
end
%%
relaPAvg(:,cycleNum:end)=[];
deltasToR1Mat(:,cycleNum:end)=[];
rPhaseMat=relaPAvg;


%%
deltasT = cycleMat(1:4,:);

deltasT=deltasT(:,~isnan(sum(deltasT,1)));


%%
plotDeltas(deltasT(1,:),deltasT(2:4,:),true)
title('time based')

R1L2P=rPhaseMat(4,:);%R1-L2
L2R3P=-(rPhaseMat(4,:)-rPhaseMat(2,:));%L2-R3
R1L1P=(rPhaseMat(3,:));%R1-L1
deltasP=[R1L2P;L2R3P;R1L1P]/(2*pi);
plotDeltas(deltasT(1,:),deltasP,true)
title('phase based')


%%

function [deltas0, deltas1, deltasTT, deltasToR1, dist, GDI, GDI_E] = getDataForPlots(sequence,refPoint,stanceDonePoint,cycleDuration,t,s)
%Speed variant
%Speed variant m-tripod May 2020
%from time shifts
mf_t = -(0.0025132*s-0.13029);
hm_t = -(0.003343*s-0.12209);
hf_t = hm_t+mf_t;

%from phase shifts
%mf_t = -(0.0037167*s-0.17575);
%hm_t = -(0.002158*s-0.063129);
%hf_t = hm_t+mf_t;

%mf_t = polyval(td.p_mf_t,speed);
%hf_t = polyval(td.p_hf_t,speed);


%extracting data just for this one cycle.
legsAndFrames=sequence(1:2,refPoint:stanceDonePoint);
%also check the leg at refPoint belongs to left or right tripod
%For delta calculations, tripod gait that starts with right
%tripod and tripod gait that starts with left tripod should be 
%treated the same way. If it starts with left tripod, switch
%that stance start times between left and right legs.
%If -1, left. If 1, right.
leftOrRight = sequence(3,refPoint);

%get frame number for each leg
R1frame = min(legsAndFrames(2,legsAndFrames(1,:)==1));
R2frame = min(legsAndFrames(2,legsAndFrames(1,:)==2));
R3frame = min(legsAndFrames(2,legsAndFrames(1,:)==3));
L1frame = min(legsAndFrames(2,legsAndFrames(1,:)==4));
L2frame = min(legsAndFrames(2,legsAndFrames(1,:)==5));
L3frame = min(legsAndFrames(2,legsAndFrames(1,:)==6));
if leftOrRight==-1
    R1temp=R1frame;
    R1frame=L1frame;
    L1frame=R1temp;
    
    R2temp=R2frame;
    R2frame=L2frame;
    L2frame=R2temp;
    
    R3temp=R3frame;
    R3frame=L3frame;
    L3frame=R3temp;
end

%We want to calculate R1-L2, L2-R3, R1-L1
%This is used for plotting Fig. 2b.
deltas0=[(t(R1frame)-t(L2frame))/cycleDuration;...
    (t(L2frame)-t(R3frame))/cycleDuration;...
    (t(R1frame)-t(L1frame))/cycleDuration];

%For comparing tripod vs tetrapod. Is there one gaussian
%deltasTT = [(t(L1frame)-t(R3frame))/cycleDuration;...
%     (t(L2frame)-t(R1frame))/cycleDuration;...
%     (t(R2frame)-t(L3frame))/cycleDuration];
deltasTT = [(t(R1frame)-t(L3frame))/cycleDuration;...
    (t(L1frame)-t(R3frame))/cycleDuration];

deltasToR1 = [(t(R2frame)-t(R1frame))/cycleDuration;...
    (t(R3frame)-t(R1frame))/cycleDuration;...
    (t(L1frame)-t(R1frame))/cycleDuration;
    (t(L2frame)-t(R1frame))/cycleDuration;
    (t(L3frame)-t(R1frame))/cycleDuration];

%Now get deltas for GDI
%stance start times (unit: cycle) into one array
ssCycle = t([R1frame R2frame R3frame L1frame L2frame L3frame])/cycleDuration;
deltas1 = getDeltas(ssCycle,'time','exp');
%Ideal tripod
ssCycleIdealTri = [0 0.5 0 0.5 0 0.5];
deltas1IdealTri = getDeltas(ssCycleIdealTri,'time','ideal');
%Get GDI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%deltas1=wrapNearIdeal(deltas1,deltas1IdealTri,'time');
[deltas1, deltas1IdealTri]=wrapExpAndIdeal(deltas1,deltas1IdealTri,'time');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GDI = getGDI(deltas1,deltas1IdealTri);

%Now do the same for tetrapods.
tetra1 = [0 2/3 1/3 1/3 0 2/3];
tetra2 = [0 2/3 1/3 2/3 1/3 0];
meta = [0 4/6 2/6 3/6 1/6 5/6];
realTri = [0 0.5+mf_t hf_t 0.5 mf_t 0.5+hf_t];
deltasIdealTetra1 = getDeltas(tetra1,'time','ideal');
deltasIdealTetra2 = getDeltas(tetra2,'time','ideal');
deltasIdealMeta = getDeltas(meta,'time','ideal');
deltasIdealRealTri = getDeltas(realTri,'time','ideal');
%deltas1T1=wrapNearIdeal(deltas1,deltasIdealTetra1,'time');
%deltas1T2=wrapNearIdeal(deltas1,deltasIdealTetra2,'time');
%deltas1M =wrapNearIdeal(deltas1,deltasIdealMeta,'time');
%deltas1RT =wrapNearIdeal(deltas1,deltasIdealRealTri,'time');
[deltas1T1, deltasIdealTetra1]=wrapExpAndIdeal(deltas1,deltasIdealTetra1,'time');
[deltas1T2, deltasIdealTetra2]=wrapExpAndIdeal(deltas1,deltasIdealTetra2,'time');
[deltas1M, deltasIdealMeta]=wrapExpAndIdeal(deltas1,deltasIdealMeta,'time');
[deltas1RT, deltasIdealRealTri]=wrapExpAndIdeal(deltas1,deltasIdealRealTri,'time');
GDI_T1 = getGDI(deltas1T1,deltasIdealTetra1);
GDI_T2 = getGDI(deltas1T2,deltasIdealTetra2);
GDI_M = getGDI(deltas1M,deltasIdealMeta);
GDI_RT = getGDI(deltas1RT,deltasIdealRealTri);
dist = [GDI; GDI_T1; GDI_T2; GDI_M; GDI_RT];
%angles = [acosd(GDI); acosd(GDI_T1); acosd(GDI_T2)];

%%%%%Calculate error%%%%%
sampling_error= mean(diff(t));
T = cycleDuration;

triE = getErrorProp(sampling_error,T,GDI,deltas1,deltas1IdealTri);
tet1E = getErrorProp(sampling_error,T,GDI_T1,deltas1T1,deltasIdealTetra1);
tet2E = getErrorProp(sampling_error,T,GDI_T2,deltas1T2,deltasIdealTetra2);
mtriE = getErrorProp(sampling_error,T,GDI_M,deltas1M,deltasIdealMeta);

GDI_E = [triE; tet1E; tet2E; mtriE];

end

function finalError = getErrorProp(sampling_error,T,d,deltaExp, deltaIdeal)
finalError = 2*sampling_error/(T*d)*sum(abs(deltaExp-deltaIdeal).*(1+abs(deltaExp)/2));
end

function plotDeltas(v,d,timeBased)
if timeBased
    ylims_a = [-0.4 0.4];
    yticks_a = [-0.25 0 0.25];
    
    ylims_b = [-0.9 -0.1];
    yticks_b = [-0.75 -0.5 -0.25];
    
    border=-0.5;
    
    %titlename='time based';
else
    ylims_a = [-0.4 0.4]*2*pi;
    yticks_a = [-0.25 0 0.25]*2*pi;
    
    ylims_b = [-0.9 -0.1]*2*pi;
    yticks_b = [-0.75 -0.5 -0.25]*2*pi;
    
    border=-0.5*2*pi;
    
    %titlename='phase based';
end

%notOutLier = ~isoutlier(d,'quartiles',2);
notOutLier = true(size(d,1),size(d,2));
 
%Remove outliers
v1 = v(notOutLier(1,:));
v2 = v(notOutLier(2,:));
v3 = v(notOutLier(3,:));
d1 = d(1,notOutLier(1,:));
d2 = d(2,notOutLier(2,:));
d3 = d(3,notOutLier(3,:));
disp(num2str(median(d1,'omitnan')))
disp(num2str(median(d2,'omitnan')))

figure
subplot(3,1,1)
hold on
scatter(v1,d1,20,209/255*[1 1 1],'filled');
plot([0 30],[0 0],'k--')
x=v1(~isnan(v1));y=d1(~isnan(d1)); doRobustFit(x,y);
%plotBoxes(v1,d1,'k')
hold off
xlim([0 28])
ylim(ylims_a)
%xticks([0 10 20 30])
yticks(yticks_a)
ylabel('R1-L2')
title(['n=' num2str(length(find(~isnan(v2+d2))))])
subplot(3,1,2)
hold on
scatter(v2,d2,20,209/255*[1 1 1],'filled');
plot([0 30],[0 0],'k--')
x=v2(~isnan(v2));y=d2(~isnan(d2)); doRobustFit(x,y);
%plotBoxes(v2,d2,'k')
hold off
xlim([0 28])
ylim(ylims_a)
%xticks([0 10 20 30])
yticks(yticks_a)
ylabel('L2-R3')
title(['n=' num2str(length(find(~isnan(v2+d2))))])
subplot(3,1,3)
hold on
scatter(v3,d3,20,209/255*[1 1 1],'filled');
plot([0 30],[1 1]*border,'k--')
%plotBoxes(v3,d3,'k')
x=v3(~isnan(v3));y=d3(~isnan(d3)); doRobustFit(x,y);
hold off
% plot(fit(tVelANorm(~isnan(tDeltaNorm(1,:)))',tDeltaNorm(1,~isnan(tDeltaNorm(1,:)))','poly2'));
% plot(fit(tVelANorm(~isnan(tDeltaNorm(2,:)))',tDeltaNorm(2,~isnan(tDeltaNorm(2,:)))','poly2'));
% plot(fit(tVelANorm(~isnan(tDeltaNorm(3,:)))',tDeltaNorm(3,~isnan(tDeltaNorm(3,:)))','poly2'));
%hold off
xlim([0 28])
ylim(ylims_b)
%xticks([0 10 20 30])
yticks(yticks_b)
ylabel('R1-L1')
title(['n=' num2str(length(find(~isnan(v3+d3))))])
end

function [y2, outliers]=doRobustFit(x,y)
[b,stat] = robustfit(x,y,'bisquare'); y2 = polyval(fliplr(b'),x);
disp(['[offest slope]: [' num2str(b') ']'])
disp(['[offest slope] p-vals: [' num2str(stat.p') ']'])

plot(x,y2,'r')

% outliers = find(abs(stat.resid)>stat.mad_s);
% disp(['% of potential outliers: ' num2str(length(outliers)/length(x)*100)])
% scatter(x(outliers),y(outliers),20,'r','filled');
end

function [tri_md, tet_md]=getMenIdx(legsVec)
numLegsInStance = sum(legsVec,1);

largeDur = size(legsVec,2);
smallDur = length(find(numLegsInStance==3));

tri_md = smallDur/largeDur;

tet_md = length(find(numLegsInStance==4))/largeDur;
end