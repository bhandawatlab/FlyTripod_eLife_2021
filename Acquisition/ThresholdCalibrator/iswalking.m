% 1) Take 3d matrix (short video)
% 2) for left, center, and right portions of the video, take average value
% of the elements in each portion and then multiply a predefined factor.
%
% @Chanwoo Chun, <cc2465@cornell>


function [m1, m2, m3]=iswalking(F)

global w rROI;

factor = 0.000001;

leftavg = sum(sum(sum(F(1:w,1:100,:))))*factor/size(F,3);
midavg = sum(sum(sum(F(1:w,rROI-50:rROI+50,:))))*factor/size(F,3);
rightavg = sum(sum(sum(F(1:w,(size(F,2)-99):size(F,2),:))))*factor/size(F,3);

m1=leftavg
m2=midavg
m3=rightavg


%if  (initsum >= threshold) && (finalsum >= threshold)
%    y = true;
%else
%    y = false;
%end


