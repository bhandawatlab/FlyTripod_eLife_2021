% This callback function updates the displayed frame and perform analysis
% on current frame.
% @Chanwoo Chun, <cc2465@cornell.edu>

function PreviewUpdate(obj,event,hImage)

global prevcount savecount tempframes diffarray vid rROI w;

%Threshold
t1 = 1.5;

numavg = 2; %Number of frames to average. Need this to smooth out noise.
%Change the value to appropriate value for your application, if needed. 

w = 1264/2-rROI(2);

% Display the current image frame. 
set(hImage, 'CData', event.Data(1:w,:));

new = event.Data;
new = new(1:w,:);

n=prevcount;
tempframes(:,:,n)= new;

if rem(n,numavg)==0
    %Quantify motion in a given range of frames
    diffarray(n-(numavg-1):n)=diffavg(tempframes(:,:,n-(numavg-1):n));
    diffarray = diffarray(1:n);
    
    %If the quantified value is above a certain threshold:
    if diffarray(n) >= t1
        
        %start acquisition
        trigger(vid);
        disp('Recording...');
        [frames, timeStamp] = getdata(vid);
        s = struct('timeStamp', timeStamp, 'frames', frames);
        framesq = squeeze(frames);
        
        %Check if the fly walked across the frame. If so, save the vid.
        if iswalking(framesq) == true            
            saver(s);
            savecount=savecount+1;
        else
            disp('Not worthy. Resume monitering... Press any key to stop.');
            clear frames;
        end
        start(vid);
    end
end

subplot(5,1,2);
plot(diffarray);

drawnow


if prevcount==1000
    prevcount=0;
    clear tempframes
    clear diffarray
end

prevcount=prevcount+1;

if savecount>=600
    stoppreview(vid);
end

