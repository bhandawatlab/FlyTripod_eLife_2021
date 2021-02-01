%This code finds initial CoM position of the fly in the sideview and the
%bottomview, and also finds a bounding boxes around the CoM points where 
%feature points will be found in another separate function.
%If a body length is 1, we say that CoM is located 0.4 away from the head
%and 0.6 away from the tail. Since we are inputting the head and tail
%locations into this function, CoM position at the bottomview (x,y) and x
%position of CoM at the sideview are easy to get. However, z position of
%CoM is not straightforward. This is explained later in the code.
%
%@Chanwoo Chun <cc2465@cornell.edu>

function [comSide,comBottom, recSide, recBottom] = findCOM(objectFrame, endFrame, background,floorY,wallUpper,wallLower,bottomHead,bottomTail)
%objectFrame is an initial image frame.
%background is a background image.
%rest of the inputs are determined by user in another function that comes
%before this function.

%background subtract image, and then separate the bottomview and sideview.
subtracted=objectFrame - background;
sideBottom = subtracted;
sideMask = zeros(size(sideBottom));
sideMask(1:floorY,:) = 1;
bottomMask = zeros(size(sideBottom));
bottomMask(wallUpper:wallLower,:)=1;
sideView = sideBottom.*uint8(sideMask);
bottomView = sideBottom.*uint8(bottomMask);

%get binarized image
%If there are more than two blobs after binarization, the code will run
%into an error. Make sure that the threshold is appropriate for imbinarize.
sideView = imbinarize(sideView,0.04);
bottomView = imbinarize(bottomView,0.04);
bottomView = bwareaopen(bottomView,100,4);

sideView = imopen(sideView,strel('disk',8));%remove legs
bottomView = imopen(bottomView,strel('disk',8));

%check if there is only one blob in each plot.
figure
subplot(2,1,1)
imshow(sideView)
subplot(2,1,2)
imshow(bottomView)
drawnow

%get centroids, length of minor axes, and orientations of the blobs.
s  = regionprops(sideView ,'centroid','MinorAxisLength','Orientation');
b  = regionprops(bottomView ,'centroid','MinorAxisLength','Orientation');

%Two blocks of for-loops shown below are used to find the best rectangles
%centered at CoM that will ensure maximum number of feature points. The for
%loops iterate through different sizes of the rectangles and see how many
%feature points can be detected inside the rectangular region.
i=30;
rectangle('Position',[s.Centroid(:,1)-4*i/2,s.Centroid(:,2),4*i,i],'EdgeColor','r');
%%sideView
maxPoints = 1;
for i=2:8:round(s.MinorAxisLength)+15

    sideRegion=[s.Centroid(:,1)-2*i/2,s.Centroid(:,2)-i/2,2*i,i];
    rectangle('Position',[s.Centroid(:,1)-2*i/2,s.Centroid(:,2)-i/2,2*i,i],'EdgeColor','w');
    try
    points = detectMinEigenFeatures(subtracted,'ROI',sideRegion);
    catch
        continue
    end
    
    if size(points.Location,1)>maxPoints
        maxPoints = size(points.Location,1);
        recSide = sideRegion;
    end

end  

i=40;
rectangle('Position',[b.Centroid(:,1)-i/2,b.Centroid(:,2)-i/2,i,i],'EdgeColor','w');
%%bottomView
maxPoints = 1;
for i=2:8:round(b.MinorAxisLength)+15

    bottomRegion=[b.Centroid(:,1)-2*i/2,b.Centroid(:,2)-2*i/2,2*i,2*i];
    rectangle('Position',[b.Centroid(:,1)-2*i/2,b.Centroid(:,2)-2*i/2,2*i,2*i],'EdgeColor','w');
    %needs to be improved -> rectangle('Position',[bottomHead(1),bottomHead(2),abs(bottomHead(1) - bottomTail(1)),abs(bottomHead(2) - bottomTail(2))],'EdgeColor','w');
    try
    points = detectMinEigenFeatures(subtracted,'ROI',bottomRegion);
    catch
    end
    
    if size(points.Location,1)>maxPoints
        maxPoints = size(points.Location,1);
        recBottom = bottomRegion;
    end

end      


%%%%%%%%%%%%%%%%%%%%%%%%%%adjust COM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comBottom=(bottomHead - bottomTail)*0.6+bottomTail;
comSide(1) = comBottom(1);
%find z position of the CoM by referring to the slope of the fly in
%sideview.
comSide(2) = tan(deg2rad(-s.Orientation))*(comBottom(1)-s.Centroid(1))+s.Centroid(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,1,1)
imshow(sideView);
hold on
plot(comSide(1),comSide(2), 'b*')
hold off
subplot(2,1,2)
imshow(bottomView)
hold on
plot(comBottom(1),comBottom(2), 'b*')
hold off
drawnow