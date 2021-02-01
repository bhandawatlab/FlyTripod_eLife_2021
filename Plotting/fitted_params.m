%This code plots the fitted parameters 
%Figure 6-S1
%@Chanwoo Chun, Jan. 31, 2021, <cc2465@cornell.edu>

pathToMatlab = strcat(extractBefore(pwd,'MATLAB'),'MATLAB');
addpath(genpath(pathToMatlab));

filepath = ['..' filesep '..' filesep 'Data' filesep 'StepData' filesep];
load([filepath 'stepData2000.mat']);
data3=stepData;

filepath = ['..' filesep '..' filesep 'Data' filesep 'StepData' filesep 'proper_fit' filesep];
load([filepath 'stepData2.mat']);
data=stepData;

g = 9807;
omitThis = zeros(1,10000);
Midx=1;
Nidx=1;
omitI=1;
ka=nan(1,length(data)); ks=ka; Ra=ka; R=ka; m=ka; L=ka; flynums=ka; rm=ka;
omega = ka; dRdt0 = ka; AP = ka; R0 = ka;
indexing = ka;
for i=1:length(data)
    %l=data3{i}.tripod.L;
    %comz=data3{i}.tripod.yMid;
    if ~isfield(data{i},'ARSLIP')
       continue
    end
    
    flynums(i) = data{i}.source.flynum;
    %L(i)=l;
    %m(i)=data{i}.source.weight/1000;
    %R(i)=data{i}.source.legLength;
    Ra(i)=data{i}.ARSLIP.Rnat;
    ka(i)=data{i}.ARSLIP.Ka;
    ks(i)=data{i}.ARSLIP.Ks;
    omega(i) = data{i}.ARSLIP.omega;
    dRdt0(i) = data{i}.ARSLIP.dRdt0;
    AP(i) = data{i}.ARSLIP.AP;
    R0(i) = data{i}.ARSLIP.R0;
    
    %rm(i)=comz;
    indexing(i)=i;
%     if  data{i}.PTR <25
%         omitThis(omitI) = i;
%         omitI = omitI+1;
%     end
end

figure
subplot(1,7,1)
verticalScatter(Ra)
title('Rnat')
subplot(1,7,2)
verticalScatter(ka)
title('ka')
subplot(1,7,3)
verticalScatter(ks)
title('ks')
subplot(1,7,4)
verticalScatter(omega)
title('omega')
subplot(1,7,5)
verticalScatter(dRdt0)
title('dRdt0')
subplot(1,7,6)
verticalScatter(AP)
title('AP')
subplot(1,7,7)
verticalScatter(R0)
title('R0')

function verticalScatter(data)
hold on
scatter(ones(size(data,2),1)*1, data, 50,'k','.', 'jitter','on', 'jitterAmount', 0.15,'LineWidth',3);
plot([0.5 1.5],[1 1]*median(data,'omitnan'),'r')
hold off
end