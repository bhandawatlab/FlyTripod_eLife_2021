%Generates gait maps
%
%@Chanwoo Chun, Aug. 21, 2019, <cc2465@cornell.edu>
%This code does the following:
% 1) Saves gait maps of all data
addpath(genpath(['..' filesep '..' filesep '..' filesep '..' filesep 'FlyLocomotion']))
shotDir= dir(['..' filesep '..' filesep '..' filesep '**' filesep 'shot.mat']);

for i = 1:length(shotDir)
    disp(num2str(i/length(shotDir)*100))
    shotName = strcat(shotDir(i).folder,filesep,shotDir(i).name);
    load(shotName);
    
    %[gender, strain, weight, legLength, flynum] = getFlyInfo(shotDir(i).folder);
    
    t = shot.timeStamp';
    
    dxdt=shot.com.bottom.dxdt*24/1984;
    dydt=shot.com.bottom.dydt*24/1984;
    comVel = (dxdt.^2+dydt.^2).^(1/2);
    
    height = shot.com.side.location(:,2)*24/1984;
    
    legs = [shot.leg.R1.';shot.leg.R2.';shot.leg.R3.';shot.leg.L1.';shot.leg.L2.';shot.leg.L3.'];   
    
    %legs = 1-legs;
    %valid=shot.validity.valid;
    %legs(:,~valid) = nan;
    
    legs(legs==0) = nan;
    
    
    legs(1,:) = legs(1,:)*10; %R1
    legs(2,:) = legs(2,:)*6; %R2
    legs(3,:) = legs(3,:)*8; %R3
    legs(4,:) = legs(4,:)*7; %L1
    legs(5,:) = legs(5,:)*9; %L2
    legs(6,:) = legs(6,:)*5; %L3
    
    lw = 10;
    figure
    ax1=subplot(3,1,1);
    hold on
    for j = 1:6
        plot(t,legs(j,:),'k','linewidth',lw)
    end
    hold off
    ylim([3 12])
    %pbaspect([3 1 1])
    
    ax2=subplot(3,1,2);
    plot(t,comVel)
    ylabel('speed (mm/s)')
    xlabel('time (sec)')
    
    ax3=subplot(3,1,3);
    plot(t,height)
    ylabel('height (mm)')
    xlabel('time (sec)')    
    linkaxes([ax1 ax2 ax3],'x')
    
    drawnow
    
    avgSpeed = mean(comVel);
    
    filePath = ['temp_data' filesep 'speed' num2str(avgSpeed) '_fileNum' num2str(i) '.png'];
    set(gcf,'visible','off')
    saveas(gcf,filePath)
    close
    
end
