%Function for realtime video analysis.
%
% @Chanwoo Chun, <cc2465@cornell>


function PreviewUpdateC(obj,event,hImage)
% This callback function updates the displayed frame and the histogram.

global prevcount fincount B C D E F vid w rROI;

nfa = 2; %Number of frames to average
%comparison ROI [240, 450, 1504, 182] failed...
%comparison ROI [480, 450, 1024, 182]
w = 1264/2-rROI(2);

% Display the current image frame. 
%set(hImage, 'CData', event.Data(yi:yf,xi:xf));
set(hImage, 'CData', event.Data(1:w,:));

new = event.Data;
new = new(1:w,:);

n=prevcount;
B(:,:,n)= new;


if rem(n,nfa)==0
    C(n-(nfa-1):n)=diffavg(B(:,:,n-(nfa-1):n));
    diffVal = C(n-(nfa-1):n)
    C = C(1:n);
    [m1, m2, m3] = iswalking(B(:,:,n-(nfa-1):n));
    
    D(n-(nfa-1):n) = m1;
    D = D(1:n);
    
    E(n-(nfa-1):n) = m2;
    E = E(1:n);
    
    F(n-(nfa-1):n) = m3;
    F = F(1:n);
    
end

subplot(5,1,2);
plot(C);

subplot(5,1,3);
plot(D);
ylim([0,1]);

subplot(5,1,4);
plot(E);
ylim([0,1]);

subplot(5,1,5);
plot(F);
ylim([0,1]);

drawnow

if prevcount==1000
    prevcount=0;
    clear B;
    clear C;
    clear D;
    clear E;
end

prevcount=prevcount+1;

if fincount>=300
    stoppreview(vid);
end

