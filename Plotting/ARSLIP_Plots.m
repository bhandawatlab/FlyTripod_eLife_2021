%This code generates:
%
% 1) Scatter plots of Ka, Ks, gammaA, gammaS, and spread (L/R) of males and
% females.
% 2) Scatter plot that shows RMSE between best fit of the models (ARSLIP
% and SLIP) and experimental data. (Figure 6C)
% 3) Plot that compares actual trends of the best fit model and
% experimental data in xy coordinate system. (Figure 6A and 6B)
%
%@Chanwoo Chun, Jan. 31, 2021, <cc2465@cornell.edu>

pathToMatlab = strcat(extractBefore(pwd,'MATLAB'),'MATLAB');
addpath(genpath(pathToMatlab));

filepath = ['..' filesep '..' filesep 'Data' filesep 'StepData' filesep];
load([filepath 'stepData2000.mat']);
data3=stepData;

filepath = ['..' filesep '..' filesep 'Data' filesep 'StepData' filesep 'proper_fit' filesep];
load([filepath 'stepData2.mat']);
data=stepData;

g = 9807;  %mm/s^2

Ka = NaN(1,length(data));
r=Ka;
Ks=Ka;
gammaA=Ka;
gammaS = Ka;
flynums = Ka;
yMid = Ka;
Lmm=Ka;
fr=Ka;
gRMSE = NaN(length(data),2);
sRMSE = NaN(length(data),2);
expRMSE = NaN(length(data),2);
Bidx = 1;
Feidx = 1;
Maidx = 1;
B = zeros(1,10000);
Fe = zeros(10000);
Ma = zeros(10000);
flynum = nan;

for i=1:length(data)
    if ~isfield(data{i},'ARSLIP') || data{i}.PTR<25
        continue
    end
    
    if strcmp(data{i}.source.strain,'Berlin K')
        B(Bidx) = i;
        Bidx=Bidx+1;
    end
    
    if strcmp(data{i}.source.gender,'female')
        Fe(Feidx) = i;
        Feidx=Feidx+1;
    else
        Ma(Maidx) = i;
        Maidx=Maidx+1;
    end
    
    flynums(i) = data{i}.source.flynum;
    
    Ka(i)=data{i}.ARSLIP.Ka;
    r(i)=data{i}.ARSLIP.Rnat;
    Ks(i)=data{i}.ARSLIP.Ks;
    m=data{i}.source.weight/1000;
    fr(i) = (data{i}.avgVel)^2/(r(i)*g);
    
    gammaA(i) = Ka(i)/(m*g*r(i));
    gammaS(i) = Ks(i)*r(i)/(m*g);
    
    %[allPos, validLegs] = getTripodPos(data{i});

    
%     [Lmm(i), yMid(i), ~] = getSpreadAndHeight(data{i});
    Lmm(i)=data3{i}.tripod.L;
    yMid(i)=data3{i}.tripod.yMid;
        
    yexp = data{i}.com;
    ysim = data{i}.ARSLIP.ysim;        
    gRMSE(i,:)=sqrt(mean((ysim-yexp).^2,1));
    ysim = data3{i}.SLIP.ysim;       
    sRMSE(i,:)=sqrt(mean((ysim-yexp).^2,1));

    % Get errors
    [errorSide, errorBottom] = getError(data{i});
    expRMSE(i,:)= sqrt(mean([errorBottom errorSide].^2,1));
    
     % Make Fit trajectory comparison figure. i=271 is a good example.
    if i == 271 %isfield(data{i},'AsySLIP')  && isfield(data{i},'GSLIP')
    validLegs = [data{i}.leg.R1';data{i}.leg.R2';data{i}.leg.R3';data{i}.leg.L1';data{i}.leg.L2';data{i}.leg.L3'];
    plotTrajComp(errorSide, errorBottom, data3{i}, i, validLegs, 'SLIP')
    plotTrajComp(errorSide, errorBottom, data{i}, i, validLegs, 'ARSLIP')
    end

 
end
%Indices for Berlin K, female, and male flies. Here, we are removing
%excessive elements in the array.
B(Bidx:end)=[];
Fe(Feidx:end)=[];
Ma(Maidx:end)=[];

normSpread = Lmm./yMid;

flyNumsBF=flynums(intersect(B,Fe));
flyNumsBM=flynums(intersect(B,Ma));
flyNumCountBF=flyCount(flyNumsBF);
flyNumCountBM=flyCount(flyNumsBM);


flyNumCount=flyCount(flynums);


figure
subplot(1,5,1)
verticalScatter(Ka,intersect(B,Fe),intersect(B,Ma))
xlabel('Berlin K')
subplot(1,5,2)
verticalScatter(Ks,intersect(B,Fe),intersect(B,Ma))
subplot(1,5,3)
verticalScatter(gammaA,intersect(B,Fe),intersect(B,Ma))
subplot(1,5,4)
verticalScatter(gammaS,intersect(B,Fe),intersect(B,Ma))
subplot(1,5,5)
verticalScatter(normSpread,intersect(B,Fe),intersect(B,Ma))


RMSEplot(sRMSE,gRMSE,271);


function verticalScatter(data,Fe,Ma)
hold on
scatter(ones(size(data(Fe),2),1)*1, data(Fe), 50,'r','.', 'jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
scatter(ones(size(data(Ma),2),1)*2, data(Ma), 50,'b','.', 'jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
plot([0.5 1.5],[1 1]*median(data(Fe),'omitnan'),'k')
plot([1.5 2.5],[1 1]*median(data(Ma),'omitnan'),'k')
hold off
xlim([0 3])
ylabel(inputname(1))
title(['p=' num2str(ranksum(data(Fe),data(Ma))) ' ma=' num2str(size(data(Ma),2)) ' fe=' num2str(size(data(Fe),2))])
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

function RMSEplot(aRMSE,bRMSE,exampleIndex)
aTemp = sum(aRMSE,2);
bTemp = sum(bRMSE,2);
aDataNum = length(find(~isnan(aTemp)));
bDataNum = length(find(~isnan(bTemp)));

sH = median(aRMSE(:,1),'omitnan');
gH = median(bRMSE(:,1),'omitnan');
sV = median(aRMSE(:,2),'omitnan');
gV = median(bRMSE(:,2),'omitnan');

allH=[aRMSE(:,1); bRMSE(:,1)];
allV=[aRMSE(:,2); bRMSE(:,2)];

figure
subplot(1,2,1)
hold on
aRMSEH = aRMSE(:,1);
%scatter(zeros(size(sRMSEH,1),1), sRMSEH, 50,'k','.','jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
scatter(zeros(size(aRMSEH,1),1), aRMSEH, 50,[1 0.227 0.22],'.','jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
%scatter(zeros(size(expRMSEH,1),1),expRMSEH,50, 'b.','jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
%scatter(0, sRMSE(doi,1),200,'r.');

scatter(0,aRMSE(exampleIndex,1),'m')

plot([-1 1],[1 1]*(1/(25.39*2)),'r')
plot([-1 1],[sH sH],'k')
hold off
%ylim([0 max(sRMSE([N M],1))*1.1])
ylim([min(allH)*0.9 max(allH)*1.1])
set(gca, 'YScale', 'log')
title(['Horizontal SLIP RMSE n= ' num2str(aDataNum)])

subplot(1,2,2)
hold on
bRMSEH = bRMSE(:,1);
%scatter(zeros(size(gRMSEH,1),1), gRMSEH, 50,'k','.', 'jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
scatter(zeros(size(bRMSEH,1),1), bRMSEH, 50,[0.38 0.443 1],'.', 'jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
%scatter(zeros(size(expRMSEH,1),1),expRMSEH,50, 'b.','jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
%scatter(0, gRMSE(doi,1),200,'r.');

scatter(0,bRMSE(exampleIndex,1),'m')

plot([-1 1],[1 1]*(1/(25.39*2)),'r')
plot([-1 1],[gH gH],'k')
hold off
%ylim([0 max(sRMSE([N M],1))*1.1])
ylim([min(allH)*0.9 max(allH)*1.1])
set(gca, 'YScale', 'log')
title(['Horizontal ARSLIP RMSE n= ' num2str(bDataNum)])
set(gcf, 'Position', [400 100 300 450], 'Visible', 'on')

figure
subplot(1,2,1)
hold on
aRMSEV = aRMSE(:,2);
%scatter(zeros(size(sRMSEV,1),1), sRMSEV, 50,'k','.', 'jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
scatter(zeros(size(aRMSEV,1),1), aRMSEV, 50,[1 0.227 0.22],'.', 'jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
%scatter(zeros(size(expRMSEV,1),1),expRMSEV,50, 'b.','jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
%scatter(0, sRMSE(doi,2),200,'r.');

scatter(0,aRMSE(exampleIndex,2),'m')

plot([-1 1],[sV sV],'k')
plot([-1 1],[1 1]*(1/(25.39*2)),'r')
hold off
%ylim([0 max(sRMSE([N M],2))*1.1])
ylim([min(allV)*0.9 max(allV)*1.1])
set(gca, 'YScale', 'log')
title(['Vertical SLIP RMSE n= ' num2str(aDataNum)])

subplot(1,2,2)
hold on
bRMSEV = bRMSE(:,2);
%scatter(zeros(size(gRMSEV,1),1), gRMSEV, 50,'k','.', 'jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
scatter(zeros(size(bRMSEV,1),1), bRMSEV, 50,[0.38 0.443 1],'.', 'jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
%scatter(zeros(size(expRMSEV,1),1),expRMSEV,50, 'b.','jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
%scatter(0, gRMSE(doi,2),200,'r.');

scatter(0,bRMSE(exampleIndex,2),'m')

plot([-1 1],[gV gV],'k')
plot([-1 1],[1 1]*(1/(25.39*2)),'r')
hold off
%ylim([0 max(sRMSE([N M],2))*1.1])
ylim([min(allV)*0.9 max(allV)*1.1])
set(gca, 'YScale', 'log')
title(['Vertical ARSLIP RMSE n= ' num2str(bDataNum)])
set(gcf, 'Position', [400 100 300 450], 'Visible', 'on')

[pHori,hHori] = ranksum(aRMSE(:,1),bRMSE(:,1))%,'Alpha',0.005)
[pVerti,hVerti] = ranksum(aRMSE(:,2),bRMSE(:,2))%,'Alpha',0.005)
end
