% This code finds the best fitting K (springy tripod's
% spring constant) to individual flies (NOT to steps). Check if leg spread
% alone can explain change in gammaA and gammaS as predicted by optimal
% fits.
%
% Figure 7C, 8
%
% @Chanwoo Chun, Jan. 31, 2021, <cc2465@cornell.edu>

pathToMatlab = strcat(extractBefore(pwd,'MATLAB'),'MATLAB');
addpath(genpath(pathToMatlab));

filepath = ['..' filesep '..' filesep 'Data' filesep 'StepData' filesep];
load([filepath 'stepData2000.mat']);
data=stepData;

filepath = ['..' filesep '..' filesep 'Data' filesep 'StepData' filesep 'proper_fit' filesep];
load([filepath 'stepData2.mat']);
datak=stepData;

g = 9807;
omitThis = zeros(1,10000);
Midx=1;
Nidx=1;
omitI=1;
ka=nan(1,length(data)); ks=ka; Ra=ka; R=ka; m=ka; L=ka; flynums=ka; rm=ka;
indexing = ka;
for i=1:length(data)
    l=data{i}.tripod.L;
    comz=data{i}.tripod.yMid;
    if isnan(comz) || isnan(l) || ~isfield(data{i},'ARSLIP')
       continue
    end
    
    flynums(i) = data{i}.source.flynum;
    L(i)=l;
    m(i)=data{i}.source.weight/1000;
    R(i)=data{i}.source.legLength;
    Ra(i)=datak{i}.ARSLIP.Rnat;
    ka(i)=datak{i}.ARSLIP.Ka;
    ks(i)=datak{i}.ARSLIP.Ks;
    rm(i)=comz;
    indexing(i)=i;
%     if  data{i}.PTR <25
%         omitThis(omitI) = i;
%         omitI = omitI+1;
%     end
end

wildFlynums = flynums;
uniqueWild = unique(flynums);
kaE = cell(1,length(uniqueWild));
ksE = cell(1,length(uniqueWild));
gaE = cell(1,length(uniqueWild));
gsE = cell(1,length(uniqueWild));
RaE = cell(1,length(uniqueWild));
IDX = cell(1,length(uniqueWild));
Kcell = cell(1,length(uniqueWild));
Rcell = cell(1,length(uniqueWild));
gammaCell = cell(1,length(uniqueWild));
LRCell = cell(1,length(uniqueWild)); 
parpool('local',4);
for j=1:length(uniqueWild)
    oneFly = wildFlynums==uniqueWild(j);
    L7=L(oneFly); ka7=ka(oneFly); ks7=ks(oneFly); Ra7=Ra(oneFly); rm7=rm(oneFly); indexing7=indexing(oneFly);
    
    R7=R(oneFly); R7=unique(R7); %this actual leg length will be used as an initial searching point.
    m7=m(oneFly); m7=unique(m7);
    
    if length(find(oneFly)) < 6
       continue 
    end
    
    L7t=L7; ka7t=ka7; ks7t=ks7; Ra7t=Ra7; rm7t = rm7;
    L7(isnan(L7t+ka7t+ks7t+Ra7t+rm7t))=[];
    ka7(isnan(L7t+ka7t+ks7t+Ra7t+rm7t))=[];
    ks7(isnan(L7t+ka7t+ks7t+Ra7t+rm7t))=[];
    Ra7(isnan(L7t+ka7t+ks7t+Ra7t+rm7t))=[];
    rm7(isnan(L7t+ka7t+ks7t+Ra7t+rm7t))=[];
    indexing7(isnan(L7t+ka7t+ks7t+Ra7t+rm7t))=[];
    
    disp([num2str(length(find(oneFly))) ' ' num2str(j/length(uniqueWild)*100)])
    
    objective = @(X) objectiveFunc(X,m7,L7,ka7,ks7,Ra7,rm7);
    k0 = 7;
    R0 = R7;
    init = [k0, R0];
    lb = [3, R0*1/5];
    ub = [16, R0*2];
    options = optimoptions('fmincon','Algorithm','interior-point','Display','iter', ...
        'GradObj','off');
    problem = createOptimProblem('fmincon','x0',init,'objective',objective,...
        'lb',lb,'ub',ub,'options',options);
    gs7 = GlobalSearch('NumTrialPoints', 200);
    [reconstructed, f] = run(gs7, problem);
    Kselected = reconstructed(1);
    Rselected = reconstructed(2);
    Kbest(j)    = Kselected;
    Rbest(j)    = Rselected;
    RMSE(j)     = f;
    massUnq(j)     = m7;
    RUnq(j)        = R7;
    
    
    ksExpect = NaN(1,length(L7));
    kaExpect = NaN(1,length(L7));
    RaExpect = NaN(1,length(L7));
    gsExpect = NaN(1,length(L7));
    gaExpect = NaN(1,length(L7));
    gs7 = NaN(1,length(L7));
    ga7 = NaN(1,length(L7));
    for i = 1:length(L7)
        [kae,kse,Rae]=solveForARSLIP(rm7(i),L7(i),Rselected,Kselected);
        ksExpect(i) = kse;
        kaExpect(i) = kae;
        RaExpect(i) = Rae;
        gsExpect(i) = kse*Rae/(m7*g);
        gaExpect(i) = kae/(m7*g*Rae);
        gs7(i) = ks7(i)*Ra7(i)/(m7*g);
        ga7(i) = ka7(i)/(m7*g*Ra7(i));
    end
    ksE{j} = [ks7;ksExpect];
    kaE{j} = [ka7;kaExpect];
    gsE{j} = [gs7;gsExpect];
    gaE{j} = [ga7;gaExpect];
    RaE{j} = [Ra7;RaExpect];
    IDX{j} = indexing7;
    Kcell{j} = Kselected*ones(1,length(ks7));
    Rcell{j} = Rselected*ones(1,length(ks7));
    
    gammaCell{j} = Kselected*Rselected/(m7*g)*ones(1,length(ks7));
    LRCell{j} = [L7; rm7];
    
    disp(num2str(Kbest))
end
kaEM=cell2mat(kaE);
ksEM=cell2mat(ksE);
gaEM=cell2mat(gaE);
gsEM=cell2mat(gsE);
RaEM=cell2mat(RaE);
IDXM=cell2mat(IDX);
KM=cell2mat(Kcell);
RM=cell2mat(Rcell);
gc_exp=cell2mat(gammaCell);
lr_exp=cell2mat(LRCell);

%%
figure
subplot(1,3,1)
hold on
scatter(gaEM(1,:),gaEM(2,:),'k','filled')
plot([min(gaEM) max(gaEM)],[min(gaEM) max(gaEM)],'r')
hold off
pbaspect([1 1 1])
Ra = corrcoef(gaEM(1,:)',gaEM(2,:)');
R2a = Ra(1,2)^2;
xlabel('gammaA ARSLIP fit')
ylabel('gammaA predicted')
title(['R2=' num2str(R2a)])

subplot(1,3,2)
hold on
scatter(gsEM(1,:),gsEM(2,:),'k','filled')
plot([min(gsEM) max(gsEM)],[min(gsEM) max(gsEM)],'r')
hold off
pbaspect([1 1 1])
Rs = corrcoef(gsEM(1,:)',gsEM(2,:)');
R2s = Rs(1,2)^2;
title(['R2=' num2str(R2s) ' n=' num2str(length(gaEM))])
xlabel('gammaS ARSLIP fit')
ylabel('gammaS predicted')
set(gcf, 'Position', [400 100 1000 500], 'Visible', 'on')

subplot(1,3,3)
hold on
scatter(RaEM(1,:),RaEM(2,:),'k','filled')
plot([min(RaEM) max(RaEM)],[min(RaEM) max(RaEM)],'r')
hold off
pbaspect([1 1 1])
Rr = corrcoef(RaEM(1,:)',RaEM(2,:)');
R2r = Rr(1,2)^2;
title(['R2=' num2str(R2r)])
xlabel('Ra ARSLIP fit')
ylabel('Ra predicted')
set(gcf, 'Position', [400 100 1000 500], 'Visible', 'on')

figure
x=(1-1./gsEM(1,:)).*RaEM(1,:);
y=(1-1./gsEM(2,:)).*RaEM(2,:);
hold on
scatter(x,y,'k','filled')
plot([min([x y]) max([x y])],[min([x y]) max([x y])],'r')
hold off
pbaspect([1 1 1])
Rrf = corrcoef(x',y');
R2rf = Rrf(1,2)^2;
title(['R2=' num2str(R2rf)])
xlabel('Rf ARSLIP fit')
ylabel('Rf predicted')
set(gcf, 'Position', [400 100 1000 500], 'Visible', 'on')


%%
minV = 0.5;
maxV = 3.5;
n=linspace(minV,maxV,100);
gc = solveGc(n);
n_exp = lr_exp(1,:)./lr_exp(2,:); 
figure
hold on
plot(n,gc,'b','linewidth',2)
%plot([minV maxV],[gamma gamma])
scatter(n_exp,gc_exp,15,'k','filled');
[f,xi] = ksdensity(n_exp); 
plot(xi,f/8+1,'r','linewidth',2);

hold off
ylabel('equ D.21')
xlabel('L/H')
title(['horizontal line is gamma. n=' num2str(length(n_exp))])

%%

delete(gcp);

function f=objectiveFunc(toBeFitted,m,L,ka,ks,Ra,rm)
g = 9807;
k=toBeFitted(1);
R=toBeFitted(2);
ksExpect = NaN(1,length(L));
kaExpect = NaN(1,length(L));
gsExpect = NaN(1,length(L));
gaExpect = NaN(1,length(L));
RaExpect = NaN(1,length(L));
gs = NaN(1,length(L));
ga = NaN(1,length(L));
parfor i = 1:length(L)
[kae,kse,Rae]=solveForARSLIP(rm(i),L(i),R,k)
ksExpect(i) = kse;
kaExpect(i) = kae;
RaExpect(i) = Rae;
gsExpect(i) = kse*Rae/(m*g);
gaExpect(i) = kae/(m*g*Rae);
gs(i) = ks(i)*Ra(i)/(m*g);
ga(i) = ka(i)/(m*g*Ra(i));
end
slopeS = mean(abs(gsExpect./gs-1));
slopeA = mean(abs(gaExpect./ga-1));
slopeR = mean(abs(RaExpect./Ra-1));
f=slopeS+slopeA+slopeR;

% Rs = corrcoef(gs',gsExpect'); R2s = Rs(1,2)^2;
% Ra = corrcoef(ga',gaExpect'); R2a = Ra(1,2)^2;
% f = -(R2s+R2a);
end

function [ka,ks,Ra]=solveForARSLIP(rm,L_mine,R,k)
% L=2*L_mine;
% 
% ka=4*k*L^2*rm^2*R/(L^2+4*rm^2)^(3/2);
% ks=k*(3-4*L^2*R/(L^2+4*rm^2)^(3/2));
% Ra=k/ks*(R-rm*(3-4*R/(L^2+4*rm^2)^(1/2)))+rm;

H=rm;
L=L_mine;
ka = k*L^2*H^2*2*R/((L^2+H^2)^(3/2));
ks = k*(3-2*L^2*R/((L^2+H^2)^(3/2)));
Ra = k/ks*(R-H*(3-2*R/((L^2+H^2)^(1/2))))+rm;
end

function gc = solveGc(n)
gc=((n.^2+1).^(3/2))./(2*n.^2);
end
