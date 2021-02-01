function plotTrajComp(errorSide, errorBottom, data, i, validLegs, model)
if isnan(validLegs)
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%Plotting the actual fit%%%%%%%%%%%%%%%%%%%%
time = data.time;
height = data.com(:,2);
distance = data.com(:,1);
velocity = data.vel;
PTR = data.PTR;

%errorVel = sqrt(movsum(errorBottom.^2,2))./gradient(time);
errorTime = std(diff(time));
termA = movsum(errorBottom.^2,2)./(gradient(distance)).^2;
termB = errorTime.^2./gradient(time).^2;
termC = abs(gradient(distance)./gradient(time));
errorVel = termC.*(termA+termB).^(1/2);

loSide = height - errorSide;
hiSide = height + errorSide;

loBottom = distance - errorBottom;
hiBottom = distance + errorBottom;

loVel = velocity - errorVel; %calculating used basic error propagation function.
hiVel = velocity + errorVel;

validLegs(validLegs == 0) = NaN;

ylimMiddle = (max(height)+min(height))/2;
ylimRange = [ylimMiddle-0.16 ylimMiddle+0.15];

%             ysim = NaN(length(time),2);
%             ysim(data{i}.pureTriStarts:data{i}.pureTriEnds,:)=data{i}.GSLIP.ysim;
if strcmp(model,'SLIP')
ysim = data.SLIP.ysim;
APpct = data.SLIP.APpct;
r = data.SLIP.Rnat;
omega = data.SLIP.omega;
gammaA = [];
gammaS = data.SLIP.Ks*r/(data.source.weight/1000*9807);
Rvel = data.SLIP.dRdt0;    
elseif strcmp(model,'ARSLIP')
ysim = data.ARSLIP.ysim;
APpct = data.ARSLIP.APpct;
r = data.ARSLIP.Rnat;
omega = data.ARSLIP.omega;
gammaA = data.ARSLIP.Ka/(data.source.weight/1000*9807*r);
gammaS = data.ARSLIP.Ks*r/(data.source.weight/1000*9807);
Rvel = data.ARSLIP.dRdt0;
elseif strcmp(model,'ST')
ysim = data.SpringyTripod2.ysim;
APpct = 0;%data.SpringyTripod2.APpct;
r = data.SpringyTripod2.Rnat;
omega = data.SpringyTripod2.omega;
gammaA = [];
gammaS = data.SpringyTripod2.K*r/(data.source.weight/1000*9807);
Rvel = data.SpringyTripod2.dRdt0;
elseif strcmp(model,'IP')
ysim = data.IP.ysim;
APpct = data.IP.APpct;
r = data.IP.R;
omega = data.IP.omega;
gammaA = [];
gammaS = [];
Rvel = [];
end

figure
subplot(3,1,1)
hpSide = patch([time; time(end:-1:1); time(1)], [loSide; hiSide(end:-1:1); loSide(1)], 'r');
hold on;
hlSide = line(time,height);
plot(time,ysim(:,2),'--k','LineWidth',1)

plot(time,validLegs(1,:)*(ylimMiddle-0.1-0.01),'k','LineWidth',2)
plot(time,validLegs(5,:)*(ylimMiddle-0.1-0.01*2),'k','LineWidth',2)
plot(time,validLegs(3,:)*(ylimMiddle-0.1-0.01*3),'k','LineWidth',2)
plot(time,validLegs(4,:)*(ylimMiddle-0.1-0.01*4),'k','LineWidth',2)
plot(time,validLegs(2,:)*(ylimMiddle-0.1-0.01*5),'k','LineWidth',2)
plot(time,validLegs(6,:)*(ylimMiddle-0.1-0.01*6),'k','LineWidth',2)

hold off;
set(hpSide, 'facecolor', [1 0.8 0.8], 'edgecolor', 'none');
set(hlSide, 'color', 'r');
ylim(ylimRange); %New Y range for qualitative trend evaluation.
%ylim([0.5 1.2]) %Original Y range.
xlabel('Time (sec)')
ylabel('Height (mm)')
title({['AP:',num2str(APpct),'%']...
    ,['Rnat:',num2str(r),'mm,Ga:',num2str(gammaA),',Gs:',num2str(gammaS)]...
    ,['Omega:',num2str(omega),'rad/s',',Rvel:',num2str(Rvel),'mm/s']})


subplot(3,1,2)
hpBottom = patch([time; time(end:-1:1); time(1)], [loBottom; hiBottom(end:-1:1); loBottom(1)], 'r');
hold on;
hlBottom = line(time,distance);
plot(time,ysim(:,1),'--k','LineWidth',1)
hold off;
set(hpBottom, 'facecolor', [1 0.8 0.8], 'edgecolor', 'none');
set(hlBottom, 'color', 'r');
xlabel('Time (sec)')
ylabel('Distance (mm)')

subplot(3,1,3)
velGSLIP =gradient(ysim(:,1))./gradient(time);
hpVel = patch([time; time(end:-1:1); time(1)], [loVel; hiVel(end:-1:1); loVel(1)], 'r');
hold on;
hlVel = line(time,velocity);
plot(time,velGSLIP,'--k','LineWidth',1)
hold off;
set(hpVel, 'facecolor', [1 0.8 0.8], 'edgecolor', 'none');
set(hlVel, 'color', 'r');
xlabel('Time (sec)')
ylabel('Speed (mm/s)')

fig = gcf;
%set(fig, 'Position', [100 100 470 632], 'Visible', 'off')
set(fig, 'Position', [100 50 470 948], 'Visible', 'on')
%         cd(strrep(data{i}.source.videoDir,' ','\'))
%         saveas(fig,strcat('Fit',num2str(data{i}.source.startFrame),'.png'))
%         cd(matdir)
shotDir=data.source.videoDir;
[filepath,vidName,~] = fileparts(shotDir);
[filepath,datName,~] = fileparts(filepath);
newName = [num2str(round(PTR)) 'pctIDX' num2str(i)];
if strcmp(model,'SLIP')
newPath = ['..' filesep '..' filesep 'Data' filesep 'FitPlot' filesep 'SLIP']; %EntireAsymGSLIP
elseif strcmp(model,'ARSLIP')
newPath = ['..' filesep '..' filesep 'Data' filesep 'FitPlot' filesep 'ARSLIP']; %EntireAsymGSLIP
elseif strcmp(model,'ST')
newPath = ['..' filesep '..' filesep 'Data' filesep 'FitPlot' filesep 'ST']; %EntireAsymGSLIP
elseif strcmp(model,'IP')
newPath = ['..' filesep '..' filesep 'Data' filesep 'FitPlot' filesep 'IP'];
end
mkdir(newPath);
saveas(fig,fullfile(newPath,newName),'png');%'png');%'epsc')
%saveas(fig,newName,'epsc');%'png');%'epsc')

end