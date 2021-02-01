%This folder looks at instantenous yaw and beta (slip) data to output
%logical indices that tells whether frames are valid (true) for analysis or
%not (false).
%
%@Chanwoo Chun, <cc2465@cornell.edu>

function [smallBeta, smallYaw, valid] = validity(shot)
%make smallBeta validity column
beta=shot.com.bottom.beta;
smallBeta = zeros(size(shot.frame)); %for now 0 means good 1 means bad
badFrameB = find(abs(beta)>10);
smallBeta(badFrameB) = 1;
smallBeta=denoise(smallBeta.',30).';
smallBeta=ones(size(smallBeta))-smallBeta; %flip so that 1 means good

%make smallYaw validity column
instYaw=shot.com.bottom.instYaw;
cn = 20;
convIYaw = conv(instYaw,ones(cn,1),'same');
convIYaw = abs(convIYaw);
smallYaw = ones(size(shot.frame));
badFrameY = find(convIYaw>8); %10<-5
smallYaw(badFrameY)=0;

%below denoises the binary array
a=smallYaw';
for i = 1:2    
    a=denoise(a,20);
    a=ones(size(a))-a;
end
smallYaw=a';

%"valid" is where fly is walking straight.
valid=shot.validity.notTouching.*smallBeta.*smallYaw;

 figure
 plot(smallBeta*9);
 hold on
 plot(smallYaw*10);
 plot(shot.validity.notTouching*11);
 plot(convIYaw);
 plot(valid);
 hold off

