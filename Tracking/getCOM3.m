% This code tracks CoM through time, using KLT tracker.

% First version: 12/6/2017
% @Chanwoo Chun <cc2465@cornell.edu>

function y=getCOM3(M,frBegin,frEnd,background)

mmPerPix = 1984/24; %mm

timeInt = M.data.timeStamp(frBegin:frEnd);%-s.data.timeStamp(frBegin);

%Video player ... Think about ROI
videoPlayer = vision.VideoPlayer('Position',[100,100,1984,255]);

%Select general tracking region, COM, floor and wall
%Advance to the frame where I can see the fly's whole body.

objectFrame = M.data.frames(:,:,:,frBegin);%now at frBegin(th) frame
endFrame =  M.data.frames(:,:,:,frEnd);
figure; imshow(objectFrame);
%%Get the side view region where you want to track.

%initial orientation
bottomHead = getPosition(impoint);
bottomTail = getPosition(impoint);

floorY = getPosition(impoint);
floorY = floorY(:,2);
wallUpper = getPosition(impoint);
wallUpper = wallUpper(:,2);
wallLower = getPosition(impoint);
wallLower = wallLower(:,2);

[comSide,comBottom,sideRegion,bottomRegion] = findCOM(objectFrame,endFrame, background,floorY,wallUpper,wallLower,bottomHead,bottomTail);
%zoneWidth = 40;

tracker = vision.PointTracker('NumPyramidLevels',3,'MaxBidirectionalError',1);

for i = 1:2
    n = 1;
    if i == 1
        %COM side view
        points = detectMinEigenFeatures(objectFrame,'ROI',sideRegion);
        comTrack = comSide;
    else
        %COM bottom view
        points = detectMinEigenFeatures(objectFrame,'ROI',bottomRegion);
        comTrack = comBottom;
    end
    
    points = points.Location;
    initialize(tracker,points,objectFrame);

    oldPoints = points;
    
    oriNetwork(:,1) = oldPoints(:,1) - comTrack(:,1);
    oriNetwork(:,2) = oldPoints(:,2) - comTrack(:,2);
    
    markers = cell(100000,1);
    markers{1} = points;
    
    errr = zeros(frBegin-frEnd);
    theta_recovered = zeros(frBegin-frEnd);
    
    %figure
    for j=frBegin+1:frEnd
          frame = M.data.frames(:,:,:,j);%-background; %now at frBegin+1

          [points, isFound] = step(tracker,frame);

          visiblePoints = points(isFound,:);
          oldInliers = oldPoints(isFound,:);
          
          if size(visiblePoints, 1) >= 2 % need at least 2 points

              % Estimate the geomatric transformation between the old points
              ...and the new points and eliminate outliers
              [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
              oldInliers, visiblePoints, 'similarity', 'MaxDistance', 1);

              % Apply the transformation to the COM!
              comTrack(n+1,:) = transformPointsForward(...
              xform, comTrack(n,:));

              %newNetwork(:,1) = visiblePoints(:,1) - comTrack(n+1,1);
              %newNetwork(:,2) = visiblePoints(:,2) - comTrack(n+1,2);
              
              %oriNetwork = oriNetwork(isFound,:);
              
              
              out = insertMarker(frame, comTrack(:,:), '+', ...
                  'Color', 'red');

              %estimating scale and rotation
              Tinv  = xform.invert.T;
              ss = Tinv(2,1);
              sc = Tinv(1,1);
              %scale_recovered = sqrt(ss*ss + sc*sc)
              theta_recovered(n+1,:) = -atan2(ss,sc)*180/pi;
              
              
              markers{n+1} = visiblePoints';
              errr=0;
              

              %Display tracked points
              out = insertMarker(out, visiblePoints, '+', ...
                  'Color', 'white');

              % Reset the points
             oldPoints = visiblePoints;
             setPoints(tracker, oldPoints);  
          else
              disp('Incomplete');    
          end

          step(videoPlayer,out);

          drawnow;
       n=n+1;
    end
    markers(n+1:end,:) = [];
    
    
    if i == 1
        comTrack(:,1) = abs(comTrack(:,1)-comTrack(1,1));
        
        %in the future, select the left and right end of the walkway
        ...and then use that information to find the height.
        com.side.location(:,2) = floorY - comTrack(:,2);
        com.side.location(:,1) = comTrack(:,1);
        com.side.error = errr;
        com.side.markers = markers;
        com.side.uncty = cumsum(abs(errr));
        com.side.instPitch = theta_recovered;
        com.side.pitch = cumsum(theta_recovered);
    else
        rawbottom = comTrack;
        
        %This is commented out since this causes trouble in future
        %analysis. This code can still be used for plotting bottom
        %position.>
        %comTrack(:,1) = abs(comTrack(:,1)-comTrack(1,1));
        
        %=wallLower - comTrack(:,2); the below code was like this. but
        %commented out for similar reason.
        com.bottom.location(:,2) = comTrack(:,2);
        com.bottom.location(:,1) = comTrack(:,1);
        com.bottom.error = errr;
        com.bottom.markers = markers;
        com.bottom.uncty = cumsum(abs(errr));
        com.bottom.instYaw = theta_recovered;
        com.bottom.yaw = cumsum(theta_recovered);
    end
    
    clear comTrack errr theta_recovered oriNetwork oriPoints markers
    release(tracker);
    
end

release(videoPlayer);

com.worldinfo.floorY = floorY;
com.worldinfo.wallUpper = wallUpper;
com.worldinfo.wallLower = wallLower;
com.worldinfo.bodyLength = pdist([bottomHead; bottomTail]);
com.worldinfo.bottomHead = bottomHead;
com.worldinfo.bottomTail = bottomTail;

%get the absolute yaw angle
headXY = bottomHead-bottomTail;
com.worldinfo.initYaw = -atan2(headXY(2),headXY(1))*180/pi; %need negative sign in the front since the (0,0) is at the upper left corner
com.bottom.yaw = com.worldinfo.initYaw+com.bottom.yaw;
%com.side.velocity = [0 diff(com.side.location(:,1))]./s.data.timeStamp(frBegin:frEnd);


dxdt=diff(rawbottom(:,1))./diff(timeInt);
dydt=diff(-rawbottom(:,2))./diff(timeInt);
dxdt=[dxdt(1); dxdt]; %we are "extrapolating" the first value by repeating
dydt=[dydt(1); dydt];
velAng = atan2(dydt,dxdt)*180/pi;

com.bottom.dxdt=dxdt;
com.bottom.dydt=dydt;

%sideslip angle arcos((v(dot)yaw_vector)/v)
%com.bottom.beta = acosd((cosd(com.bottom.yaw).*dxdt+sind(com.bottom.yaw).*dydt)./((dxdt.^2+dydt.^2).^(1/2)));
com.bottom.beta = velAng-com.bottom.yaw;

figure
      subplot(4,1,1);
      plot(com.side.location(:,2));
      subplot(4,1,2);
      plot(com.side.location(:,1));
      subplot(4,1,3);
      plot(com.bottom.yaw);
      subplot(4,1,4);
      plot(com.bottom.beta);
      drawnow

      figure
      subplot(2,1,1)
      plot(cumsum(com.side.error))
      subplot(2,1,2)
      plot(cumsum(com.bottom.error))      
      
y = com;
    end
    