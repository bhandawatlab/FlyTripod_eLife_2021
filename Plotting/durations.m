%This code generates the swing and stance duration vs speed plot.
%Figure 1D
%
%@Chanwoo Chun, Jan. 31, 2021, <cc2465@cornell.edu>

addpath(genpath(['..' filesep '..' filesep '..' filesep 'FlyLocomotion']))
shotDir= dir(['..' filesep '..' filesep '**' filesep 'shot.mat']);

stC=cell(1,10000);
swC=cell(1,10000);
cnt=1;
for i = 1:length(shotDir)
    [gender, strain, weight, legLength, flynum] = getFlyInfo(shotDir(i).folder);
    
    load([shotDir(i).folder filesep shotDir(i).name]);
    
    valid=shot.validity.valid;
    validLegs = [shot.leg.R1.';shot.leg.R2.';shot.leg.R3.';shot.leg.L1.';shot.leg.L2.';shot.leg.L3.'];
    
    %The below four lines were needed because of ONE strange data we had.
    %In the data, one leg was tracked outside the valid portion of the
    %data. This messed up the analysis. Therefore, the data is omitted.
    summed=sum(validLegs,1);
    if any(summed(logical(1-valid)))
        continue
    end
    
    time = shot.timeStamp';%*1000;
    dxdt=shot.com.bottom.dxdt*24/1984;
    dydt=shot.com.bottom.dydt*24/1984;
    comVel = (dxdt.^2+dydt.^2).^(1/2);
    
    %determin starts and ends of the valid portion
    [validStarts, validEnds] = getStartsAndEnds(valid);
    for j=1:length(validStarts)
        vStart=validStarts(j);
        vEnd=validEnds(j);
        
        %Loop through the legs.
        for h=1:6
            comVelPortion=comVel(vStart:vEnd);
            leg=validLegs(h,vStart:vEnd);
            
            %Get stance durations of a given leg in a given valid portion.
            [legStarts, legEnds] = getStartsAndEnds(leg);
            [stStarts, stEnds]=removeInvalidStartsAndEnds(legStarts,legEnds,1,size(leg,2),'option2');
            stDurs=time(stEnds)-time(stStarts);
            
            stSpeed=nan(1,length(stDurs));
            for k = 1:length(stDurs)
                stSpeed(k) = mean(comVelPortion(stStarts(k):stEnds(k)));
            end
            stC{cnt}=[stSpeed;stDurs;flynum*ones(1,length(stDurs))];
            
            %Get swing durations of a given leg in a given valid portion.
            [invLegStarts, invLegEnds] = getStartsAndEnds(1-leg);
            [swStarts, swEnds]=removeInvalidStartsAndEnds(invLegStarts,invLegEnds,1,size(leg,2),'option2');
            swDurs=time(swEnds)-time(swStarts);
            
            swSpeed=nan(1,length(swDurs));
            for k = 1:length(swDurs)
                swSpeed(k) = mean(comVelPortion(swStarts(k):swEnds(k)));
            end
            swC{cnt}=[swSpeed;swDurs;flynum*ones(1,length(swDurs))];
            
            cnt=cnt+1;
        end
    end
end
stC(:,cnt:end) = [];
swC(:,cnt:end) = [];
st = cell2mat(stC);
sw = cell2mat(swC);

st=st(:,sum(st,1)>0);
sw=sw(:,sum(sw,1)>0);

%Getting stance data count for each fly.
flynumsSt = st(3,:);
flyNumCountSt=flyCount(flynumsSt);

%Getting swing data count for each fly.
flynumsSw = sw(3,:);
flyNumCountSw=flyCount(flynumsSw);

%%

ft1 = fittype('a/x+b');
f1=fit(st(1,:)',st(2,:)',ft1)
f2=fit(sw(1,:)',sw(2,:)','poly1')
fitlm(sw(1,:),sw(2,:))
figure
hold on
scatter(st(1,:),st(2,:),10,[255, 204, 121]/255,'filled');
scatter(sw(1,:),sw(2,:),10,[122, 201, 255]/255,'filled');
ax=gca;
ax.YLim = [0 .25];
h1=plot(f1,'r');
h2=plot(f2,'b');
set(h1,'Color',[255, 157, 0]/255,'LineWidth',3)
set(h2,'Color',[0, 151, 255]/255,'LineWidth',3)
set(gca,'TickDir','out')
hold off
ylabel('time (s)')
xlabel('speed (mm/s)')
