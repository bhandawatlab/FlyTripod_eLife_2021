function y = getLegs2(shot, binVid, legTips, analysisdir, addedFirst, addedLast)

validExcptBetaYaw = shot.validity.notTouching;
valid = shot.validity.valid;
%validExcptYaw = shot.validity.notTouching.*shot.validity.smallBeta;
scbLocation = shot.com.bottom.location;
scbYaw = shot.com.bottom.yaw;

%ExtrapolationExtrapolationExtrapolationExtrapolationExtrapolationExtrapolationExtrapolationExtrapolationExtrapolationExtrapolation
%valid = [valid(1)*ones(addedFirst,1);valid;valid(end)*ones(addedLast,1)];
%valid = [valid(1)*ones(addedFirst,1);validExcptYaw;valid(end)*ones(addedLast,1)];

x = cumsum(ones(length(scbLocation(:,1)),1))+addedFirst;
xq = cumsum(ones(addedFirst,1));
startExtrapolated = interp1(x,scbLocation(:,1),xq,'linear','extrap');
xq = cumsum(ones(addedLast,1))+x(end);
endExtrapolated = interp1(x,scbLocation(:,1),xq,'linear','extrap');
scbLocationTemp(:,1) = [startExtrapolated;scbLocation(:,1);endExtrapolated];

x = cumsum(ones(length(scbLocation(:,2)),1))+addedFirst;
xq = cumsum(ones(addedFirst,1));
startExtrapolated = interp1(x,scbLocation(:,2),xq,'linear','extrap');
xq = cumsum(ones(addedLast,1))+x(end);
endExtrapolated = interp1(x,scbLocation(:,2),xq,'linear','extrap');
scbLocationTemp(:,2) = [startExtrapolated;scbLocation(:,2);endExtrapolated];

 scbLocation = scbLocationTemp;
% scbLocation
figure
plot(scbLocation(:,1))
figure
plot(scbLocation(:,2))


x = cumsum(ones(length(scbYaw),1))+addedFirst;
xq = cumsum(ones(addedFirst,1));
startExtrapolated = interp1(x,scbYaw,xq,'linear','extrap');
xq = cumsum(ones(addedLast,1))+x(end);
endExtrapolated = interp1(x,scbYaw,xq,'linear','extrap');
scbYaw = [startExtrapolated;scbYaw;endExtrapolated];
%scbYaw = [scbYaw(1)*ones(addedFirst,1);scbYaw;scbYaw(end)*ones(addedLast,1)];
%ExtrapolationExtrapolationExtrapolationExtrapolationExtrapolationExtrapolationExtrapolationExtrapolationExtrapolationExtrapolation


tdValid(1,1,:)=valid;
bodyLength = shot.com.worldinfo.bodyLength;

binVid = squeeze(binVid);
legTipsV = squeeze(legTips);

%%perform row by row multiplication.
%legTipsV=bsxfun(@times,legTipsV,tdValid);

[legLabel, num] = bwlabeln(legTipsV);

stanceLocation = cell(1,num);
allPositions=zeros(size(binVid,3),2,num);

tframes=size(binVid,3);

leg.R1 = zeros(tframes,1);
leg.R2 = zeros(tframes,1);
leg.R3 = zeros(tframes,1);
leg.L1 = zeros(tframes,1);
leg.L2 = zeros(tframes,1);
leg.L3 = zeros(tframes,1);
leg.rawPositions.R1 = zeros(tframes,2);
leg.rawPositions.R2 = zeros(tframes,2);
leg.rawPositions.R3 = zeros(tframes,2);
leg.rawPositions.L1 = zeros(tframes,2);
leg.rawPositions.L2 = zeros(tframes,2);
leg.rawPositions.L3 = zeros(tframes,2);
leg.legPositionPlot.R1 = zeros(tframes,2);
leg.legPositionPlot.R2 = zeros(tframes,2);
leg.legPositionPlot.R3 = zeros(tframes,2);
leg.legPositionPlot.L1 = zeros(tframes,2);
leg.legPositionPlot.L2 = zeros(tframes,2);
leg.legPositionPlot.L3 = zeros(tframes,2);

figure

label=11;
for label = 1:num
    
    %get logical (binary) video where only the leg of interest is 1.
    legRaw = double(legLabel == label);
    legPresent = squeeze(sum(sum(legRaw,1),2)~=0);
    startToEnd = find(legPresent)';
    legPositionRaw=zeros(size(legRaw,3),2);
    legPositionPlot=legPositionRaw;
    comPosition = scbLocation; %%%
    comPosition(:,2) = -(comPosition(:,2)-shot.com.worldinfo.wallUpper+1);
    comPosition = comPosition.*[legPresent,legPresent];
    
    %comPositionPlot = [comPosition(:,2),comPosition(:,1)];
    bodyAngle = scbYaw.*[legPresent,legPresent];
    %figure
    for f = startToEnd
        s = regionprops(legRaw(:,:,f),'centroid');
        
        legPositionRaw(f,:)=s.Centroid;
        legPositionPlot(f,:) = legPositionRaw(f,:);
        legPositionPlot(f,2) = -legPositionPlot(f,2);
        %change coordinate (flip)
        %legPositionPlot(f,:) = [legPositionRaw(f,2),legPositionRaw(f,1)];
        
        %further change (translation)
        legPositionPlot(f,:) = legPositionPlot(f,:)-comPosition(f,:);
        
        %final change (rotate)
        theta = deg2rad(-(bodyAngle(f)-90));
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        legPositionPlot(f,:) = (R*legPositionPlot(f,:).').';
        
        %       imshow(mat2gray(legRaw(:,:,f)))
        %       drawnow
    end
    
    forwardMove = legPositionPlot(:,2);
    forwardMove = mean(forwardMove(forwardMove~=0));
    lateralMove = legPositionPlot(:,1);
    lateralMove = mean(lateralMove(lateralMove~=0));
    distance = (forwardMove^2+lateralMove^2)^(1/2);
    di = 0;
    rivalEraser = ones(size(legRaw,3),1); %By default, we do not want to erase.
    
    if lateralMove < 0 %On the plot, left. For the fly, RIGHT.
        if forwardMove >= 60 % R1
            
            if any(leg.R1(legPresent))==true
                rivalLabel = max(leg.R1(legPresent)); %find who is taking the spot.
                rivalEraser = squeeze(sum(sum(double(legLabel == rivalLabel),1),2)==0); %Prepare eraser
                
                fM = leg.legPositionPlot.R1(legPresent,2);
                fM = mean(fM(fM~=0));
                lM = leg.legPositionPlot.R1(legPresent,1);
                lM = mean(lM(lM~=0));
                di = (fM^2+lM^2)^(1/2);
            end
            
            if  distance > di
                %delete the rival footprint and add the new one.
                leg.R1 = leg.R1.*rivalEraser + double(legPresent)*label;
                leg.rawPositions.R1 = leg.rawPositions.R1.*[rivalEraser rivalEraser] + legPositionRaw;
                leg.legPositionPlot.R1 = leg.legPositionPlot.R1.*[rivalEraser rivalEraser] + legPositionPlot;
            end
        elseif forwardMove <= -60 %R3
            if any(leg.R3(legPresent))==true
                rivalLabel = max(leg.R3(legPresent)); %find who is taking the spot.
                rivalEraser = squeeze(sum(sum(double(legLabel == rivalLabel),1),2)==0); %Prepare eraser
                fM = leg.legPositionPlot.R3(legPresent,2);
                fM = mean(fM(fM~=0));
                lM = leg.legPositionPlot.R3(legPresent,1);
                lM = mean(lM(lM~=0));
                di = (fM^2+lM^2)^(1/2);
            end
            
            if  distance > di
            leg.R3 = leg.R3.*rivalEraser + double(legPresent)*label;
            leg.rawPositions.R3  = leg.rawPositions.R3.*[rivalEraser rivalEraser]  + legPositionRaw;
            leg.legPositionPlot.R3  = leg.legPositionPlot.R3.*[rivalEraser rivalEraser]  + legPositionPlot;
            end
        else %R2
            if any(leg.R2(legPresent))==true
                rivalLabel = max(leg.R2(legPresent)); %find who is taking the spot.
                rivalEraser = squeeze(sum(sum(double(legLabel == rivalLabel),1),2)==0); %Prepare eraser
                fM = leg.legPositionPlot.R2(legPresent,2);
                fM = mean(fM(fM~=0));
                lM = leg.legPositionPlot.R2(legPresent,1);
                lM = mean(lM(lM~=0));
                di = (fM^2+lM^2)^(1/2);
            end
            
            if  distance > di
            leg.R2 = leg.R2.*rivalEraser + double(legPresent)*label;
            leg.rawPositions.R2  = leg.rawPositions.R2.*[rivalEraser rivalEraser]  + legPositionRaw;
            leg.legPositionPlot.R2  = leg.legPositionPlot.R2.*[rivalEraser rivalEraser]  + legPositionPlot;
            end
        end
    else %On the plot, right. For the fly, LEFT.
        if forwardMove >= 60 % L1
            if any(leg.L1(legPresent))==true
                rivalLabel = max(leg.L1(legPresent)); %find who is taking the spot.
                rivalEraser = squeeze(sum(sum(double(legLabel == rivalLabel),1),2)==0); %Prepare eraser
                fM = leg.legPositionPlot.L1(legPresent,2);
                fM = mean(fM(fM~=0));
                lM = leg.legPositionPlot.L1(legPresent,1);
                lM = mean(lM(lM~=0));
                di = (fM^2+lM^2)^(1/2);
            end
            
            if  distance > di
            leg.L1 = leg.L1.*rivalEraser + double(legPresent)*label;
            leg.rawPositions.L1 = leg.rawPositions.L1.*[rivalEraser rivalEraser] + legPositionRaw;
            leg.legPositionPlot.L1 = leg.legPositionPlot.L1.*[rivalEraser rivalEraser] + legPositionPlot;
            end
        elseif forwardMove <= -60 %L3
            if any(leg.L3(legPresent))==true
                rivalLabel = max(leg.L3(legPresent)); %find who is taking the spot.
                rivalEraser = squeeze(sum(sum(double(legLabel == rivalLabel),1),2)==0); %Prepare eraser
                fM = leg.legPositionPlot.L3(legPresent,2);
                fM = mean(fM(fM~=0));
                lM = leg.legPositionPlot.L3(legPresent,1);
                lM = mean(lM(lM~=0));
                di = (fM^2+lM^2)^(1/2);
            end
            
            if  distance > di
            leg.L3 = leg.L3.*rivalEraser + double(legPresent)*label;
            leg.rawPositions.L3  = leg.rawPositions.L3.*[rivalEraser rivalEraser]  + legPositionRaw;
            leg.legPositionPlot.L3  = leg.legPositionPlot.L3.*[rivalEraser rivalEraser]  + legPositionPlot;
            end
        else %L2
            if any(leg.L2(legPresent))==true
                rivalLabel = max(leg.L2(legPresent)); %find who is taking the spot.
                rivalEraser = squeeze(sum(sum(double(legLabel == rivalLabel),1),2)==0); %Prepare eraser
                fM = leg.legPositionPlot.L2(legPresent,2);
                fM = mean(fM(fM~=0));
                lM = leg.legPositionPlot.L2(legPresent,1);
                lM = mean(lM(lM~=0));
                di = (fM^2+lM^2)^(1/2);
            end
            
            if  distance > di
            leg.L2 = leg.L2.*rivalEraser + double(legPresent)*label;
            leg.rawPositions.L2  = leg.rawPositions.L2.*[rivalEraser rivalEraser]  + legPositionRaw;
            leg.legPositionPlot.L2  = leg.legPositionPlot.L2.*[rivalEraser rivalEraser]  + legPositionPlot;
            end
        end
    end
    
    allPositions(:,:,label) = legPositionPlot;
    
    plot(allPositions(:,1,label)./bodyLength,allPositions(:,2,label)./bodyLength,lateralMove./bodyLength,forwardMove./bodyLength,'o')
    hold on
    drawnow
        
end

leg.R1=double(leg.R1~=0);
leg.R2=double(leg.R2~=0);
leg.R3=double(leg.R3~=0);
leg.L1=double(leg.L1~=0);
leg.L2=double(leg.L2~=0);
leg.L3=double(leg.L3~=0);

%UNEXTRAPOLATEUNEXTRAPOLATEUNEXTRAPOLATEUNEXTRAPOLATEUNEXTRAPOLATEUNEXTRAPOLATEUNEXTRAPOLATE
binVid=binVid(:,:,1+addedFirst:end-addedLast);
legTipsV=legTipsV(:,:,1+addedFirst:end-addedLast);
legLabel=legLabel(:,:,1+addedFirst:end-addedLast);

leg.R1=leg.R1(1+addedFirst:end-addedLast);
leg.R2=leg.R2(1+addedFirst:end-addedLast);
leg.R3=leg.R3(1+addedFirst:end-addedLast);
leg.L1=leg.L1(1+addedFirst:end-addedLast);
leg.L2=leg.L2(1+addedFirst:end-addedLast);
leg.L3=leg.L3(1+addedFirst:end-addedLast);


leg.rawPositions.R1 = leg.rawPositions.R1(1+addedFirst:end-addedLast,:);
leg.rawPositions.R2 = leg.rawPositions.R2(1+addedFirst:end-addedLast,:);
leg.rawPositions.R3 = leg.rawPositions.R3(1+addedFirst:end-addedLast,:);
leg.rawPositions.L1 = leg.rawPositions.L1(1+addedFirst:end-addedLast,:);
leg.rawPositions.L2 = leg.rawPositions.L2(1+addedFirst:end-addedLast,:);
leg.rawPositions.L3 = leg.rawPositions.L3(1+addedFirst:end-addedLast,:);



leg.legPositionPlot.R1 = leg.legPositionPlot.R1(1+addedFirst:end-addedLast,:);
leg.legPositionPlot.R2 = leg.legPositionPlot.R2(1+addedFirst:end-addedLast,:);
leg.legPositionPlot.R3 = leg.legPositionPlot.R3(1+addedFirst:end-addedLast,:);
leg.legPositionPlot.L1 = leg.legPositionPlot.L1(1+addedFirst:end-addedLast,:);
leg.legPositionPlot.L2 = leg.legPositionPlot.L2(1+addedFirst:end-addedLast,:);
leg.legPositionPlot.L3 = leg.legPositionPlot.L3(1+addedFirst:end-addedLast,:);

%UNEXTRAPOLATEUNEXTRAPOLATEUNEXTRAPOLATEUNEXTRAPOLATEUNEXTRAPOLATEUNEXTRAPOLATEUNEXTRAPOLATE



%VALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATE
%perform row by row multiplication.
%legTipsV=bsxfun(@times,legTipsV,tdValid);
leg.R1=bsxfun(@times,leg.R1,validExcptBetaYaw);
leg.R2=bsxfun(@times,leg.R2,validExcptBetaYaw);
leg.R3=bsxfun(@times,leg.R3,validExcptBetaYaw);
leg.L1=bsxfun(@times,leg.L1,validExcptBetaYaw);
leg.L2=bsxfun(@times,leg.L2,validExcptBetaYaw);
leg.L3=bsxfun(@times,leg.L3,validExcptBetaYaw);


leg.rawPositions.R1=bsxfun(@times,leg.rawPositions.R1,validExcptBetaYaw);
leg.rawPositions.R2=bsxfun(@times,leg.rawPositions.R2,validExcptBetaYaw);
leg.rawPositions.R3=bsxfun(@times,leg.rawPositions.R3,validExcptBetaYaw);
leg.rawPositions.L1=bsxfun(@times,leg.rawPositions.L1,validExcptBetaYaw);
leg.rawPositions.L2=bsxfun(@times,leg.rawPositions.L2,validExcptBetaYaw);
leg.rawPositions.L3=bsxfun(@times,leg.rawPositions.L3,validExcptBetaYaw);



leg.legPositionPlot.R1=bsxfun(@times,leg.legPositionPlot.R1,validExcptBetaYaw);
leg.legPositionPlot.R2=bsxfun(@times,leg.legPositionPlot.R2,validExcptBetaYaw);
leg.legPositionPlot.R3=bsxfun(@times,leg.legPositionPlot.R3,validExcptBetaYaw);
leg.legPositionPlot.L1=bsxfun(@times,leg.legPositionPlot.L1,validExcptBetaYaw);
leg.legPositionPlot.L2=bsxfun(@times,leg.legPositionPlot.L2,validExcptBetaYaw);
leg.legPositionPlot.L3=bsxfun(@times,leg.legPositionPlot.L3,validExcptBetaYaw);

%VALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATE


combinedVid=uint8(zeros(size(legTipsV,1),size(legTipsV,2),3,length(find(validExcptBetaYaw))));

cVidIdx = 1;
figure;
for f=find(validExcptBetaYaw)'
    
    C = imfuse(mat2gray(legTipsV(:,:,f)),mat2gray(binVid(:,:,f)),'falsecolor','Scaling','joint','ColorChannels',[2 1 2]);

    if leg.R1(f) ~= 0
    C = insertText(C,leg.rawPositions.R1(f,:),'R1','FontSize',18,'BoxColor',...
    'yellow','BoxOpacity',0.4,'TextColor','white');
    end
    
    if leg.R2(f) ~= 0
    C = insertText(C,leg.rawPositions.R2(f,:),'R2','FontSize',18,'BoxColor',...
    'yellow','BoxOpacity',0.4,'TextColor','white');
    end
    
    if leg.R3(f) ~= 0
    C = insertText(C,leg.rawPositions.R3(f,:),'R3','FontSize',18,'BoxColor',...
    'yellow','BoxOpacity',0.4,'TextColor','white');
    end
    
    
    if leg.L1(f) ~= 0
    C = insertText(C,leg.rawPositions.L1(f,:),'L1','FontSize',18,'BoxColor',...
    'yellow','BoxOpacity',0.4,'TextColor','white');
    end
    
    if leg.L2(f) ~= 0
    C = insertText(C,leg.rawPositions.L2(f,:),'L2','FontSize',18,'BoxColor',...
    'yellow','BoxOpacity',0.4,'TextColor','white');
    end
    
    if leg.L3(f) ~= 0
    C = insertText(C,leg.rawPositions.L3(f,:),'L3','FontSize',18,'BoxColor',...
    'yellow','BoxOpacity',0.4,'TextColor','white');
    end
    
    combinedVid(:,:,:,cVidIdx) = C; 
    cVidIdx = cVidIdx + 1;
end

totalTime = shot.timeStamp-shot.timeStamp(1);
validLegs = [leg.R1.';leg.R2.';leg.R3.';leg.L1.';leg.L2.';leg.L3.'];
validLegs(validLegs==0)=NaN;
fig = figure;
plot(totalTime,validLegs(1,:)*0.6,'k')
hold on
plot(totalTime,validLegs(5,:)*0.5,'k')
plot(totalTime,validLegs(3,:)*0.4,'k')
plot(totalTime,validLegs(4,:)*0.3,'k')
plot(totalTime,validLegs(2,:)*0.2,'k')
plot(totalTime,validLegs(6,:)*0.1,'k')
ylim([0 0.7])
set(findall(gca, 'Type', 'Line'),'LineWidth',10);
hold off
set(fig, 'Position', [200 200 634 171])

fixes=legEditor(combinedVid);

for fx = 1:length(fixes) 
    labelnum = legLabel(fixes{fx}.position(2),fixes{fx}.position(1),fixes{fx}.frame);
    legRaw = double(legLabel == labelnum);
    legPresent = squeeze(sum(sum(legRaw,1),2)~=0);
    startToEnd = find(legPresent)';
    legPositionRaw=zeros(size(legRaw,3),2);
    legPositionPlot=legPositionRaw;
    comPosition = shot.com.bottom.location; %%%
    comPosition(:,2) = -(comPosition(:,2)-shot.com.worldinfo.wallUpper+1);
    comPosition = comPosition.*[legPresent,legPresent];
    
    %comPositionPlot = [comPosition(:,2),comPosition(:,1)];
    bodyAngle = shot.com.bottom.yaw.*[legPresent,legPresent];
    %figure
    for f = startToEnd
        s = regionprops(legRaw(:,:,f),'centroid');
        
        legPositionRaw(f,:)=s.Centroid;
        legPositionPlot(f,:) = legPositionRaw(f,:);
        legPositionPlot(f,2) = -legPositionPlot(f,2);
        %change coordinate (flip)
        %legPositionPlot(f,:) = [legPositionRaw(f,2),legPositionRaw(f,1)];
        
        %further change (translation)
        legPositionPlot(f,:) = legPositionPlot(f,:)-comPosition(f,:);
        
        %final change (rotate)
        theta = deg2rad(-(bodyAngle(f)-90));
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        legPositionPlot(f,:) = (R*legPositionPlot(f,:).').';
    end
    
    fields = fieldnames(leg);
    for i = 6
        rawPosTemp = leg.rawPositions.(fields{i});
        rawPosTemp(rawPosTemp == 0) = 100000;
        indicator = rawPosTemp - legPositionRaw;
        if any(indicator == 0)
            indices = find(~indicator);
            L = bwlabel(leg.(fields{i}));
            erasingMask = 1-(L == L(indices(1)));
            leg.(fields{i}) = leg.(fields{i}).*erasingMask;
            leg.rawPositions.(fields{i})=bsxfun(@times,leg.rawPositions.(fields{i}),erasingMask);
            leg.legPositionPlot.(fields{i})=bsxfun(@times,leg.legPositionPlot.(fields{i}),erasingMask);
        end
            
    end
    
    indicator = leg.(fields{fixes{fx}.leg}) + double(legPresent);
    L = bwlabel(leg.(fields{fixes{fx}.leg}));
    if any(indicator==2)
        indices = find(indicator == 2);
        erasingMask = 1-(L == L(indices(1)));
        leg.(fields{fixes{fx}.leg})=leg.(fields{fixes{fx}.leg}).*erasingMask;
        leg.rawPositions.(fields{fixes{fx}.leg})=bsxfun(@times,leg.rawPositions.(fields{fixes{fx}.leg}),erasingMask);
        leg.legPositionPlot.(fields{fixes{fx}.leg})=bsxfun(@times,leg.legPositionPlot.(fields{fixes{fx}.leg}),erasingMask);
    end
    
    switch fixes{fx}.leg
        case 1
            leg.R1 = leg.R1 + double(legPresent);
            leg.rawPositions.R1 = leg.rawPositions.R1 + legPositionRaw;
            leg.legPositionPlot.R1 = leg.legPositionPlot.R1 + legPositionPlot;
        case 2
            leg.R2 = leg.R2 + double(legPresent);
            leg.rawPositions.R2 = leg.rawPositions.R2 + legPositionRaw;
            leg.legPositionPlot.R2 = leg.legPositionPlot.R2 + legPositionPlot;            
        case 3
            leg.R3 = leg.R3 + double(legPresent);
            leg.rawPositions.R3 = leg.rawPositions.R3 + legPositionRaw;
            leg.legPositionPlot.R3 = leg.legPositionPlot.R3 + legPositionPlot;            
        case 4
            leg.L1 = leg.L1 + double(legPresent);
            leg.rawPositions.L1 = leg.rawPositions.L1 + legPositionRaw;
            leg.legPositionPlot.L1 = leg.legPositionPlot.L1 + legPositionPlot;            
        case 5
            leg.L2 = leg.L2 + double(legPresent);
            leg.rawPositions.L2 = leg.rawPositions.L2 + legPositionRaw;
            leg.legPositionPlot.L2 = leg.legPositionPlot.L2 + legPositionPlot;            
        case 6
            leg.L3 = leg.L3 + double(legPresent);
            leg.rawPositions.L3 = leg.rawPositions.L3 + legPositionRaw;
            leg.legPositionPlot.L3 = leg.legPositionPlot.L3 + legPositionPlot;            
    end
          
end

combinedVid=uint8(zeros(size(legTipsV,1),size(legTipsV,2),3,length(find(validExcptBetaYaw))));

validLegs = [leg.R1.';leg.R2.';leg.R3.';leg.L1.';leg.L2.';leg.L3.'];
validLegs(validLegs==0)=NaN;
fig = figure;
plot(totalTime,validLegs(1,:)*0.6,'k')
hold on
plot(totalTime,validLegs(5,:)*0.5,'k')
plot(totalTime,validLegs(3,:)*0.4,'k')
plot(totalTime,validLegs(4,:)*0.3,'k')
plot(totalTime,validLegs(2,:)*0.2,'k')
plot(totalTime,validLegs(6,:)*0.1,'k')
ylim([0 0.7])
set(findall(gca, 'Type', 'Line'),'LineWidth',10);
hold off
set(fig, 'Position', [200 200 634 171])

cVidIdx = 1;
figure;
for f=find(validExcptBetaYaw)'
    
    C = imfuse(mat2gray(legTipsV(:,:,f)),mat2gray(binVid(:,:,f)),'falsecolor','Scaling','joint','ColorChannels',[2 1 2]);

    if leg.R1(f) ~= 0
    C = insertText(C,leg.rawPositions.R1(f,:),'R1','FontSize',18,'BoxColor',...
    'yellow','BoxOpacity',0.4,'TextColor','white');
    end
    
    if leg.R2(f) ~= 0
    C = insertText(C,leg.rawPositions.R2(f,:),'R2','FontSize',18,'BoxColor',...
    'yellow','BoxOpacity',0.4,'TextColor','white');
    end
    
    if leg.R3(f) ~= 0
    C = insertText(C,leg.rawPositions.R3(f,:),'R3','FontSize',18,'BoxColor',...
    'yellow','BoxOpacity',0.4,'TextColor','white');
    end
    
    
    if leg.L1(f) ~= 0
    C = insertText(C,leg.rawPositions.L1(f,:),'L1','FontSize',18,'BoxColor',...
    'yellow','BoxOpacity',0.4,'TextColor','white');
    end
    
    if leg.L2(f) ~= 0
    C = insertText(C,leg.rawPositions.L2(f,:),'L2','FontSize',18,'BoxColor',...
    'yellow','BoxOpacity',0.4,'TextColor','white');
    end
    
    if leg.L3(f) ~= 0
    C = insertText(C,leg.rawPositions.L3(f,:),'L3','FontSize',18,'BoxColor',...
    'yellow','BoxOpacity',0.4,'TextColor','white');
    end

    combinedVid(:,:,:,cVidIdx) = C; 
    cVidIdx = cVidIdx + 1;
end



diskLogger = VideoWriter([analysisdir filesep num2str(shot.dataNum)], 'Motion JPEG AVI');
open(diskLogger);
for i = 1:size(combinedVid,4)
    imshow(combinedVid(:,:,:,i))
    writeVideo(diskLogger,combinedVid(:,:,:,i));
end
close
close(diskLogger);



%%%Now move the analyzed data to bigYawBeta
leg.allYawBeta = leg;
%VALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATE
%perform row by row multiplication.
%legTipsV=bsxfun(@times,legTipsV,tdValid);
leg.R1=bsxfun(@times,leg.R1,valid);
leg.R2=bsxfun(@times,leg.R2,valid);
leg.R3=bsxfun(@times,leg.R3,valid);
leg.L1=bsxfun(@times,leg.L1,valid);
leg.L2=bsxfun(@times,leg.L2,valid);
leg.L3=bsxfun(@times,leg.L3,valid);


leg.rawPositions.R1=bsxfun(@times,leg.rawPositions.R1,valid);
leg.rawPositions.R2=bsxfun(@times,leg.rawPositions.R2,valid);
leg.rawPositions.R3=bsxfun(@times,leg.rawPositions.R3,valid);
leg.rawPositions.L1=bsxfun(@times,leg.rawPositions.L1,valid);
leg.rawPositions.L2=bsxfun(@times,leg.rawPositions.L2,valid);
leg.rawPositions.L3=bsxfun(@times,leg.rawPositions.L3,valid);



leg.legPositionPlot.R1=bsxfun(@times,leg.legPositionPlot.R1,valid);
leg.legPositionPlot.R2=bsxfun(@times,leg.legPositionPlot.R2,valid);
leg.legPositionPlot.R3=bsxfun(@times,leg.legPositionPlot.R3,valid);
leg.legPositionPlot.L1=bsxfun(@times,leg.legPositionPlot.L1,valid);
leg.legPositionPlot.L2=bsxfun(@times,leg.legPositionPlot.L2,valid);
leg.legPositionPlot.L3=bsxfun(@times,leg.legPositionPlot.L3,valid);

%VALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATEVALIDATE


y = leg;
end




