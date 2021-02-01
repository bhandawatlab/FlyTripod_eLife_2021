%This code performs image processing on the bottom view.
%This is used for a preperation for tracking leg movement in getLegs2
%(called by the parent function "mainTrack").
%@Chanwoo Chun <cc2465@cornell.edu>

function [binVid, legTips, notTouching, addedFirst, addedLast] = bottomViewAnalysis(M,frBegin,frEnd,wallUpper,wallLower,background)

%Here, we are adding more frames to analyse before and after the defined
%starting frame (frBegin) and ending frame (frEnd). This will be fixed back
%by getLegs2 function. We need to add the additional frames at the front
%and at the end, since there is going to be a moving window operation 
%across the frames. By appending more data at the beginning and at the
%end, the moving window operation can be operated at those portions.
addedFirst=20;%1
addedLast=20;%1
%Make sure the code is not appending excessive amount of frames.
if frBegin<=addedFirst
    addedFirst=frBegin-1;
end
frBegin=frBegin-addedFirst;

if size(M.data.frames,4)<=frEnd+addedLast
    addedLast=size(M.data.frames,4)-frEnd;
end
frEnd=frEnd+addedLast;

wallUpper = int16(wallUpper);
wallLower = int16(wallLower);

tframes = frEnd-frBegin+1;

M.data.frames = M.data.frames(:,:,:,frBegin:frEnd);
frames = squeeze(M.data.frames);

%Declaring the video matrices.
A = zeros(wallLower-wallUpper+1,size(M.data.frames, 2),size(M.data.frames, 3),tframes);
S = zeros(wallLower-wallUpper+1,size(M.data.frames, 2),size(M.data.frames, 3),tframes);
J = zeros(wallLower-wallUpper+1,size(M.data.frames, 2),size(M.data.frames, 3),tframes);
K = zeros(wallLower-wallUpper+1,size(M.data.frames, 2),size(M.data.frames, 3),tframes);

SE = strel('disk',2);
SE2 = strel('disk',6);

%Perform frame-wise image processes. At the end, J will be a video that
%only shows the tips of the legs as white blobs. K will be a video of a
%skeleton. A will be a background subtracted raw video.
figure
for i = 1:tframes
    A(:,:,i) = M.data.frames(wallUpper:wallLower,:,i) - background(wallUpper:wallLower,:);
    
    %binarize image
    S(:,:,i) = imbinarize(A(:,:,i),4); %use 4 or 7 (7 for bright vid).
    
    %remove noise
    J(:,:,i) = bwareaopen(S(:,:,i), 20);
    
    imshow(mat2gray((J(:,:,i))))
    drawnow
    
    %dilate image
    J(:,:,i)=imdilate(J(:,:,i),strel('disk',2));
    
    %get skeleton of the fly by thinning
    K(:,:,i)=bwmorph(J(:,:,i),'thin',Inf);
  
    %get endpoints of the skeletons
   J(:,:,i) = bwmorph(K(:,:,i),'endpoints');
    
    %dilate the endpoints
   J(:,:,i)=imdilate(J(:,:,i),SE2);
end

%Determine if the fly touched the wall by using the skeleton video.
%This will be used for determining a given frame is valid for future
%analysis.
K=squeeze(sum(K,2));
touchTop = find(prod(K(1:50,:),1));    
touchBottom = find(prod(K(end-50:end,:),1));

touch = union(touchTop,touchBottom);

notTouching = ones(size(K,2),1);
notTouching(touch) = 0;

a=denoise(notTouching.',50);
    
notTouching = a.';
notTouching = notTouching(1+addedFirst:end-addedLast,:);

%Before this, J showed a leg tips at all frames. However, we only need to
%see leg tips during stance phases (footprint).
%As a reminder, J was a black video where leg tips were shown as white
%blobs (objects). We know from experiments that a fly's stance phase lasts
%about 13 frames. Therefore, we are going to look at each pixel and see how
%its value changes through time. In the time series, we are only going to
%keep a group of points that are composed of more than 13 white points
%(value 1). Smaller groups will be erased (changed to zero 0). In other
%words, it is similar to denoising. After doing this for every individual
%pixel, the resulting video will only show stance phases. I wrote a
%vectorized code for denoising. This allows us to perform this
%operation fairly quickly. Remember, this is a moving window operation and
%this is why we needed to append more frames at the beginning and end the
%end.
for h = 1:size(J,1)
   for w = 1:size(J,2)
    J(h,w,:)=denoise(squeeze(J(h,w,:))',13);
   end
   h/size(J,1)*100
end
close
binVid = A;%(:,:,:,1+addedFirst:end-addedLast);
legTips = J;%(:,:,:,1+addedFirst:end-addedLast);
