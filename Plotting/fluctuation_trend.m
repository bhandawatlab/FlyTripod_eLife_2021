% This code generates CoM characteristics vs. speed
% Figure 4 and Figure 2-S2 
%
% PTR is pure tripod ratio and equals to single stance duration / tripod
% duration.
% 
%@Chanwoo Chun, Jan. 31, 2021, <cc2465@cornell.edu>

pathToMatlab = strcat(extractBefore(pwd,'MATLAB'),'MATLAB');
addpath(genpath(pathToMatlab));
filepath = ['..' filesep '..' filesep 'Data' filesep 'StepData' filesep];
load([filepath 'stepData2000.mat']);
data=stepData;

g = 9807;  %mm/s^2

nanArray = NaN(1,length(data));
PTR = nanArray;
hde = nanArray;
vde = nanArray;
hdeAbs = nanArray; 
vdeAbs = nanArray;
avgVel = nanArray;
distTravel = nanArray; %May 24, 2020
allAvgVel = nanArray; %May 24, 2020
durations = nanArray; %June 16, 2020
initVel = nanArray; %June 16, 2020
flynums = nanArray;
expRMSE = NaN(length(data),2);
expRMSEB = NaN(length(data),1);
expSpeedRMSE = NaN(length(data),1);
expSpeedRMSEB = NaN(length(data),1);

omitThis = zeros(1,10000);
heightSeries = NaN(length(data),201);
velSeries = heightSeries;
omitI = 1;
for i=1:length(data)
    time = data{i}.time;
    avgVel(i) = data{i}.avgVel;
    theTime=data{i}.theTime;
    theHeight = data{i}.height_interp/data{i}.height_interp(1)*100-100;
    hde(i) = squeeze(getHDE(data{i}.height'));
    theVel = data{i}.vel_interp/data{i}.vel_interp(1)*100-100;
    vde(i) = squeeze(getHDE(data{i}.vel'));
    PTR(i) = data{i}.PTR;
    heightSeries(i,:) = theHeight-(max(theHeight)+min(theHeight))/2;
    velSeries(i,:) = theVel-(max(theVel)+min(theVel))/2;
    
    distTravel(i) = data{i}.com(end,1); %May 24, 2020
    allAvgVel(i) = data{i}.avgVel; %May 24, 2020
    durations(i) = time(end); %June 16, 2020
    initVel(i) = data{i}.vel(1);
    
    % Get errors
    hdeAbs(i) = squeeze(getAbsHDE(data{i}.height'));
    vdeAbs(i) = squeeze(getAbsHDE(data{i}.vel'));
    bodyLength=data{i}.source.bodyLength;
    yexp = data{i}.com;
    [errorSide, errorBottom] = getError(data{i});
    expRMSE(i,:)= sqrt(mean([errorBottom errorSide].^2,1));
    errorTime = std(diff(time));
    termA = movsum(errorBottom.^2,2)./(gradient(yexp(:,1))).^2;
    termB = errorTime.^2./gradient(time).^2;
    termC = abs(gradient(yexp(:,1))./gradient(time));
    errorSpeed = termC.*(termA+termB).^(1/2);
    expSpeedRMSE(i)= sqrt(mean(errorSpeed.^2,1));
    expRMSEB(i)=mean(errorSide)*1984/24*1/bodyLength;
    expSpeedRMSEB(i)=expSpeedRMSE(i)*1984/24*1/bodyLength;
    
    if hde(i)<-0.06 %below -0.06 is probably an artifact.
        omitThis(omitI) = i;
        omitI = omitI+1;
    else
        flynums(i) = data{i}.source.flynum;
    end
end
omitThis(omitI:end)=[];
hde(omitThis)=[];
vde(omitThis)=[];
hdeAbs(omitThis)=[];
vdeAbs(omitThis)=[];
avgVel(omitThis)=[];
expRMSEB(omitThis)=[];
expSpeedRMSEB(omitThis)=[];

allAvgVel(omitThis)=[];% May 24, 2020
distTravel(omitThis)=[];% May 24, 2020
durations(omitThis)=[];% June 16,2020
initVel(omitThis)=[];% June 16,2020
%%
figure
subplot(2,1,1)
plotHrGpBox(avgVel,hde,vde,'Hr')
title('normal')
subplot(2,1,2)
plotHrGpBox(avgVel,hde,vde,'Gp')
%%
figure
plotAvgTrend(heightSeries,velSeries,PTR,'h')
figure
plotAvgTrend(heightSeries,velSeries,PTR,'v')

%%
flyNumCount=flyCount(flynums);
%% May 24, 2020

%f1=fit(allAvgVel',distTravel','poly1')

figure
hold on
scatter(allAvgVel,distTravel,20,'k','filled')
ax=gca;
%ax.YLim = [0 .25];
%h1=plot(f1,'r');
[y2, outliers] = doRobustFit(allAvgVel,distTravel);
%set(h1,'Color',[255, 157, 0]/255,'LineWidth',3)
hold off
ylabel('dist. traveled (mm)')
xlabel('avg. speed (mm/s)')
title(['N=' num2str(length(find(~isnan(distTravel))))])

%% June 16
%ft1 = fittype('a/x+b');
%ft1 = fittype('a/(x^b)');
ft1 = fittype('a*x^b+c');
f1=fit(allAvgVel',durations'*1000,ft1,'Lower',[-inf,-inf,-inf],'Upper',[inf,0,inf])

figure
hold on
scatter(allAvgVel,durations*1000,20,'k','filled')
ax=gca;
ax.YLim = [0 .1*1000];
h1=plot(f1,'r');
set(h1,'Color',[255, 157, 0]/255,'LineWidth',3)
set(gca,'TickDir','out')
hold off
ylabel('step duration (ms)')
xlabel('avg. speed (mm/s)')
%%
ft2 = fittype('x^a');
f2=fit(allAvgVel',initVel',ft2)
figure
hold on
scatter(allAvgVel,initVel,20,'k','filled')
plot([0,1]*max([allAvgVel initVel]),[0 1]*max([allAvgVel initVel]),'r','linewidth',2)
ax2=gca;
ax2.YLim = [0 30];
h2=plot(f2,'r');
set(h2,'Color',[255, 157, 0]/255,'LineWidth',3)
set(gca,'TickDir','out')
hold off
ylabel('init. speed (mm/s)')
xlabel('avg. speed (mm/s)')
daspect([1 1 1])


%%

heratio = abs(hdeAbs)./expRMSEB';
seratio = abs(vdeAbs)./expSpeedRMSEB';

figure
subplot(1,2,1)
hold on
scatter(zeros(size(heratio,2),1), heratio, 50,'k','.','jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
%boxplot(heratio)
plot([-0.5 0.5],[1 1],'k')
plot([-0.5 0.5],[1 1]*mean(heratio,'omitnan'),'r')
%plot([1 1],[0 100],'k')
%plot([mean(heratio) mean(heratio)],[0 100],'r')
hold off
title(['Hc/Exp error. ' num2str(length(find(heratio>1))) ':' num2str(num2str(length(find(heratio<1))))])
ylim([0.5 5*10^2])
set(gca, 'YScale', 'log')
subplot(1,2,2)
hold on
scatter(zeros(size(seratio,2),1), seratio, 50,'k','.','jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
%histogram(seratio)
plot([-0.5 0.5],[1 1],'k')
plot([-0.5 0.5],[1 1]*mean(seratio,'omitnan'),'r')
%plot([1 1],[0 100],'k')
%plot([mean(seratio) mean(seratio)],[0 100],'r')
hold off
title(['Sc/Exp error. ' num2str(length(find(seratio>1))) ':' num2str(num2str(length(find(seratio<1))))])
ylim([0.5 5*10^2])
set(gca, 'YScale', 'log')


%%
function [y2, outliers]=doRobustFit(x,y)
[b,stat] = robustfit(x,y,'bisquare');
y2 = polyval(fliplr(b'),x);
disp(['[offest slope]: [' num2str(b') ']'])
disp(['[offest slope] p-vals: [' num2str(stat.p') ']'])

plot(x,y2,'Color',[255, 157, 0]/255,'LineWidth',3)

outliers = find(abs(stat.resid)>stat.mad_s);
%disp(['% of potential outliers: ' num2str(length(outliers)/length(x)*100)])
%scatter(x(outliers),y(outliers),20,'r','filled');
end

%%
function plotAvgTrend(heightSeries,velSeries,PTR,deName)

if strcmp(deName,'h')
   seriesSet=heightSeries;
   ylab = 'height fluc. (mm)';
elseif strcmp(deName,'v')
   seriesSet=velSeries;
   ylab = 'vel fluc. (mm/s)';
end

seriesSet(PTR<25,:)=[];

x = linspace(1,100,201)';
y = mean(seriesSet,1)';
e = std(seriesSet,0,1)';

lo = y - e;
hi = y + e;

%Making the mimimum value 0.
lo=lo-y(1);
hi=hi-y(1);
y=y-y(1);

hp = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], 'r');
hold on
hl = line(x,y);
hold off
ylabel(ylab)
xticks([0 25 50 75 100])
title(flag)
set(hp, 'facecolor', [190 190 190]/255, 'edgecolor', 'none');
set(hl, 'color', 'k');

if strcmp(deName,'h')
    ylim([-4 7])
    %ylim([-0.01 0.02])
    %yticks([-0.01 0 0.01 0.02])
elseif strcmp(deName,'v')
    ylim([-40 100])
    %ylim([-2 6])
    %yticks([-2 -1 0 1 2 3 4 5 6])
end

end

function plotHrGpBox(speed,hde,vde,deName)
if strcmp(deName,'Hr')
   de=hde;
   ylab = 'height fluc. amp. (mm)';
elseif strcmp(deName,'Gp')
   de=vde;
   ylab = 'vel fluc. amp. (mm/s)';
end

preX = linspace(min(speed),max(speed),100);
poly2Cons = 'a*x';%'a*x^2';
fe1=fit(speed',de',poly2Cons,'StartPoint',[0.03]);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mdl = fitlm(speed,de,'Intercept',false)
edges = [0 7.5 10 12.5 15 17.5 20 22.5 30];
Y = discretize(speed,edges);
deS = de;
boxplotM=NaN(length(deS),length(edges)-1);
p=NaN(1,length(edges)-1);
for i = 1:length(edges)-1
    boxplotM(1:length(deS(Y==i)),i)=deS(Y==i)';
    %[~,p(i)] = ttest(deS(Y==i)');

    [p(i),~,~] = signrank(deS(Y==i)',0,'tail','right');
end
p
hold on
scatter(speed,de,10,'k','filled')
plot(fe1,'r')%,'MarkerSize',30,'LineWidth',1.25)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(preX,0*preX,'k','Linewidth',1)
hold off
ylim([min(de)*1.1 max(de)*1.1])
%xlim([0 max(speed)*1.1])
xlim([0 30])
ylabel(ylab)
xlabel('speed (mm/s)')
legend('off')
if strcmp(deName,'Hr')
   yticks([-0.04 0 0.04 0.08 0.12])
elseif strcmp(deName,'Gp')
   yticks([-5 0 5 10 15 20 25 30 35])
end
xticks(edges);
xticklabels({'0', '7.5', '10', '12.5', '15', '17.5', '20', '22.5', '30'})
end

function hde = getHDE(theHeight)
line = linspace(theHeight(1),theHeight(end),length(theHeight));
deHeight = theHeight-line;

maxi = max(deHeight);
mini = min(deHeight);
hde = maxi+mini;
end

function hde = getAbsHDE(theHeight)
line = linspace(theHeight(1),theHeight(end),length(theHeight));
deHeight = theHeight-line;

maxi = max(deHeight);
mini = min(deHeight);
%hde = maxi+mini;
hde = abs(maxi)+abs(mini);
end

function [errorSide, errorBottom] = getError(data)

yourComputer = extractBefore(pwd,'Dropbox (Duke Bio_Ea)');
shotDir = data.source.videoDir; %vidDir is not actually video dir...
shotDir=strrep(shotDir,'C:\Users\cc583\',yourComputer);

if ismac
    shotDir=strrep(shotDir,'C:\Users\cc583\','/Users/badooki/');
    shotDir=strrep(shotDir,'\','/');
end
try
    load([shotDir filesep 'shot.mat']);
catch
        errorSide = NaN;
        errorBottom = NaN;
        return
end

startFrame = data.source.startFrame;
endFrame = data.source.endFrame;

if isfield(shot.com.side,'uncty')
        errorSide = NaN;
        errorBottom = NaN;
        return
end
errorSide = shot.com.side.error(startFrame:endFrame)*24/1984;
errorBottom = shot.com.bottom.error(startFrame:endFrame)*24/1984;
end
