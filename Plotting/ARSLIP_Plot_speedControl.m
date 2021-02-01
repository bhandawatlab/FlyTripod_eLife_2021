%This code generates plots that shows relationship between nondimensional
%speed and optimal parameter values of ARSLIP. Also, it shows relationship
%between nondimensional speed and spread (L/R). For each fly, a regression
%line was generated individually.
%
%Figure 6D; 7D,E
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
Bidx = 1;
Feidx = 1;
Maidx = 1;
B = zeros(1,10000);
Fe = zeros(10000);
Ma = zeros(10000);

intI = zeros(1,10000);
intrsc = 1;

%May 25,2020
midOmega=Ka;
iniOmega=Ka;
midXs=Ka;
iniXs=Ka;

%Loop through all the steps.
for i=1:length(data)
    if ~isfield(data{i},'ARSLIP') || data{i}.PTR<25 || isnan(data3{i}.tripod.L)
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
    
    %Fly number is unique to each fly.
    flynums(i) = data{i}.source.flynum;
    
    %Optimized parameter values
    Ka(i)=data{i}.ARSLIP.Ka;
    r(i)=data{i}.ARSLIP.Rnat;
    Ks(i)=data{i}.ARSLIP.Ks;
    m=data{i}.source.weight/1000;
    fr(i) = (data{i}.avgVel)^2/(r(i)*g);
    
    gammaA(i) = Ka(i)/(m*g*r(i));
    gammaS(i) = Ks(i)*r(i)/(m*g);
    
    Lmm(i)=data3{i}.tripod.L;
    yMid(i)=data3{i}.tripod.yMid;
    
    %May 25,2020
    %mid=round(size(data{i}.ARSLIP.rawSolution,1));
    [~,mid]=min(abs(data{i}.ARSLIP.ysim(:,1)-data{i}.ARSLIP.AP));
    disp(mid)
    midOmega(i)=data{i}.ARSLIP.rawSolution(mid,2);
    iniOmega(i)=data{i}.ARSLIP.rawSolution(1,2);
    xs=gradient(data{i}.ARSLIP.ysim(:,1));
    midXs(i)=xs(mid);
    iniXs(i)=xs(1);
end
flyNumCount=flyCount(flynums);

normSpread = Lmm./yMid;

figure
subplot(1,3,1)
plotFlyBasis2(fr,gammaA,flynums)
xlabel('fr')
ylabel('gamma_a')
subplot(1,3,2)
plotFlyBasis2(fr,gammaS,flynums)
xlabel('fr')
ylabel('gamma_s')
subplot(1,3,3)
plotFlyBasis2(fr,normSpread,flynums)
xlabel('fr')
ylabel('L/r0')

%May 25,2020
figure
hold on
scatter(fr,midOmega,10,'r','filled')
scatter(fr,iniOmega,10,'b','filled')
hold off
legend('mid omega','init omega')
ylabel('omega')
xlabel('fr')

figure
hold on
scatter(fr,midOmega-iniOmega,10,'k','filled')
yline(0)
hold off
ylabel('mid omega - init. omega')
xlabel('fr')
figure
subplot(1,2,1)
hold on
scatter(fr,midXs-iniXs,10,'k','filled')
yline(0)
hold off
ylabel('mid x speed - init. x speed')
xlabel('fr')
subplot(1,2,2)
hold on
boxplot(midXs-iniXs)
yline(0)
hold off
ylabel('mid x speed - init. x speed')

%%%Supplementary%%%
uniqF= unique(flynums);
uniqF(isnan(uniqF))=[];
subi=1;
for i = unique(uniqF)
Li=Lmm(flynums==i);
fri=fr(flynums==i);
[frs,sortFr]=sort(fri);
Ls=Li(sortFr);
slide = 0.1;
zx=0;
subplot(5,2,subi)
hold on
for g = 1:length(frs)
    x = zx + g*slide;
    handl = plot([x x],[-Ls(g)/2 Ls(g)/2]);
    set(handl, {'color'}, {[242 110 37]/255},'linewidth',2);
end
hold off
ylim([-2 2])
xlim([0 10])
ylabel('(mm)')
subi=subi+1;
set(gcf,'renderer','painters');
end


figure
subplot(1,2,1)
plotFlyBasis2(fr,Lmm,flynums)
xlabel('fr')
ylabel('L (mm)')
subplot(1,2,2)
plotFlyBasis2(fr,yMid,flynums)
xlabel('fr')
ylabel('Midstance height (mm)')

function plotFlyBasis2(x,y,flynums)
uniqueFlyNums = unique(flynums);
thres = 0.045;
if isrow(x)
    x=x';
end

if isrow(y)
    y=y';
end
% 
% if nargout ~= 0
% r2 = NaN(length(uniqueFlyNums),1);
% p = NaN(length(uniqueFlyNums),1);
% dataN = NaN(length(uniqueFlyNums),1);
% end
% 
% colorOrder=jet(length(uniqueFlyNums));
sc = scatter(x,y);
ylimVal = ylim;
xlimVal = xlim;
set(sc,'visible','off');
idx=1;
hold on
for i=uniqueFlyNums

    fliesOfInterest=find(i==flynums);
    xoi = x(fliesOfInterest,:);
    yoi = y(fliesOfInterest,:);
    
    numberOfData = length(fliesOfInterest);
    if numberOfData == 0
        %if it is not a normal (wildtype) fly.
        continue
    end
    %figure
    %scatter(xoi,yoi,10,'filled')
    disp(['fb i=' num2str(i) ' length=' num2str(length(xoi))])
    if  numberOfData<6
        i
       continue 
    end
    
    [fit2, rsqr, pval, n]=fittingFunc(xoi,yoi,'a*x+b');
    
    
    if pval>thres
        colorChoice = [200 205 255]/255;
        thickness = 1;
        lineS = '--';
    else
        colorChoice = [128 141 255]/255;
        thickness = 2;
        lineS = '-';
    end
    
    
    %scatter(xoi,yoi,20,colorChoice,'filled')
    
    h1=plot(fit2,'r');
    %set(h1,'Color',colorChoice,'LineWidth',1)    
    set(h1,'Color',[1 1 1]*200/255,'LineWidth',thickness,'LineStyle',lineS) 
    clear hl
    idx=idx+1;
end

% foiF=find(6==flynums);
% scatter(x(foiF),y(foiF),20,'r')
foiM=find(7==flynums);
scatter(x(foiM,:),y(foiM,:),20,'k','filled')
% [fitF, ~, ~, ~]=fittingFunc(x(foiF),y(foiF),'a*x+b');
% h1=plot(fitF,'r');
[fitM, ~, pval, ~]=fittingFunc(x(foiM,:),y(foiM,:),'a*x+b');
if pval>thres
    colorChoice = [255 199 199]/255;
    thickness = 1;
    lineS = '--';
else
    colorChoice = [255 128 128]/255;
    thickness = 2;
    lineS = '-';
end
h2=plot(fitM,'k');
set(h2,'Color','k','LineWidth',thickness,'LineStyle',lineS) 
hold off
legend(gca,'off');
ylabel(inputname(2))
xlabel(inputname(1))
xlim(xlimVal)
ylim(ylimVal)
title(['n=' num2str(length(foiM))])
end


function [fit2, rsqr, pval, n]=fittingFunc(x,y,equation)
ft = fittype(equation);
if isrow(x)
    x=x';
end

if isrow(y)
    y=y';
end

xtemp=x;
ytemp=y;

x(isnan(xtemp+ytemp))=[];
y(isnan(xtemp+ytemp))=[];

n=length(x);
    
tbl = table(x,y,'VariableNames',{'x','y'});
lm = fitlm(tbl,'y~x');

pval = lm.Coefficients.pValue(2);
rsqr = lm.Rsquared.Ordinary;


fit2 = fit(x,y,ft,'Robust','on');
end

