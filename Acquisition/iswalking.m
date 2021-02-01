%Determine if a fly was walking across a certain distance in a given video.
%If fly walked across the left and middle sections, and/or right and middle
%sections, then we say fly walked significant amount.
%@Chanwoo Chun, <cc2465@cornell.edu>

function y=iswalking(F)

global w rROI;

%Threshold
threshold = 0.2;

factor = 0.000001;

%Take a sum of pixel intensity values at the left portion of the frames
initsum = sum(sum(F(1:w,1:100,:)))*factor;

%Take a sum of pixel intensity values at the middle portion of the frames
midsum = sum(sum(F(1:w,rROI-50:rROI+50,:)))*factor;

%Take a sum of pixel intensity values at the right portion of the frames
finalsum = sum(sum(F(1:w,(size(F,2)-99):size(F,2),:)))*factor;

maxini = max(initsum);
maxmid = max(midsum);
maxfinal = max(finalsum);

disp(strcat('maxini=',num2str(maxini)));
disp(strcat('maxmid=',num2str(maxmid)));
disp(strcat('maxfinal=',num2str(maxfinal)));

subplot(5,1,3)
plot(squeeze(initsum))
subplot(5,1,4)
plot(squeeze(midsum))
subplot(5,1,5)
plot(squeeze(finalsum))
drawnow

if  (any(initsum>threshold) && any(midsum>threshold))||(any(finalsum>threshold) && any(midsum>threshold))
    y = true;
else
    y = false;
end