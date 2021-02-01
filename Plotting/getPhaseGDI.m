%This code generates leg phases of "ideal" tripod using getIdealPhases
%function. Then, it calculates GDI between the experimental phase and ideal
%phases using calculateGDI function.
%
%Chanwoo Chun Jan. 20, 2019 <cc2465@cornell.edu>

function [relativeExp, phaseGDI, dist] = getPhaseGDI(shot,analyticalSignals,t,startFrame,endFrame,refLeg,speed)
%Take a portion of the entire experimental analytical signal.
%This is the time domain of the cycle of interest.
expPhases=analyticalSignals(:,startFrame:endFrame);
%When analyticalSignal is empty vector, it is when the leg position data
%was never acquired (probably because the video was too short to have
%enough data).
%"any(isnan(sum(expPhases,1)))" is there to disable getting phaseGDI when
%phases at some frames have NaN value. However, in most cases, there are
%only two ~ five frames of NaNs per cycle which is about ~40 frames. So
%having NaN does not compromise our analysis. Therefore, it is recommended
%to just comment it out. When this omission was enabled, it omitted about
%half of our data, which is not preferred.
if isempty(analyticalSignals) %|| any(isnan(sum(expPhases,1)))
    phaseGDI = NaN;
    return
end

legPos=shot.leg.legPositionPlot;

switch refLeg
    case 1
        refPos = legPos.R1(:,2);
    case 2
        refPos = legPos.R2(:,2);
    case 3
        refPos = legPos.R3(:,2);
    case 4
        refPos = legPos.L1(:,2);
    case 5
        refPos = legPos.L2(:,2);
    case 6
        refPos = legPos.L3(:,2);
end
%The parent code assigned "refFrame", which is a stance start time of the 
%ref. leg, to "startFrame" input argument. so "startFrame" is a stance start
%time. The parent code also assigned "endFrame" which is the next stance
%start time of the ref. leg to "endFrame" input argument. so "endFrame" is
%essentially a swing end time.
%So now we only need to know stance end time.
[~, ends] = getStartsAndEnds(refPos);
stanceEnd = ends(ends>startFrame & ends<endFrame);
%just in case start are more than one value (shouldn't happen though)...
if length(stanceEnd)>1
    disp('Warning: check getPhaseGDI.')
    stanceEnd = min(stanceEnd);
end
stanceDuration = t(stanceEnd)-t(startFrame);
swingDuration = t(endFrame)-t(stanceEnd);
stanceLength = abs(refPos(stanceEnd)-refPos(startFrame));

expTime = t(startFrame:endFrame)-t(startFrame);
[idealPhases, ~, idealTime] = getIdealPhases(swingDuration, stanceDuration, stanceLength, speed, 'tripod');
[relativeExp, patternScoreT] = calculateGDI(expPhases,idealPhases,expTime,idealTime); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseGDI=mean(patternScoreT,'omitnan');

[idealPhases, ~, idealTime] = getIdealPhases(swingDuration, stanceDuration, stanceLength, speed, 'tetrapod1');
[~, patternScoreTT1] = calculateGDI(expPhases,idealPhases,expTime,idealTime); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseGDI_tetrapod1=mean(patternScoreTT1,'omitnan');

[idealPhases, ~, idealTime] = getIdealPhases(swingDuration, stanceDuration, stanceLength, speed, 'tetrapod2');
[~, patternScoreTT2] = calculateGDI(expPhases,idealPhases,expTime,idealTime); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseGDI_tetrapod2=mean(patternScoreTT2,'omitnan');

%this one is outdated despite its name (m_tripod). It's not m_tripod, it's
%actually a wave1 gait
[idealPhases, ~, idealTime] = getIdealPhases(swingDuration, stanceDuration, stanceLength, speed, 'wave1');
[~, patternScore] = calculateGDI(expPhases,idealPhases,expTime,idealTime); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseGDI_m_tripod=mean(patternScore,'omitnan');

%This is m-tripod.
[idealPhases, ~, idealTime] = getIdealPhases(swingDuration, stanceDuration, stanceLength, speed, 'realTripod');
[~, patternScoreM] = calculateGDI(expPhases,idealPhases,expTime,idealTime); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseGDI_r_tripod=mean(patternScoreM,'omitnan');

dist = [phaseGDI; phaseGDI_tetrapod1; phaseGDI_tetrapod2; phaseGDI_m_tripod; phaseGDI_r_tripod];


end

function [idealPhases, idealPositions, idealTime] = getIdealPhases(swingDuration, stanceDuration, stanceLength, s, gaitType)
%td=load('tripod_delays.mat');
%mf_p = polyval(td.p_mf_p,speed);
%hf_p = polyval(td.p_hf_p,speed);

%from phase shifts
mf_p = -(0.0037167*s-0.17575);
hm_p = -(0.002158*s-0.063129);
hf_p = hm_p+mf_p;

% mf_p = 0.083;
% hf_p = 0.083;

%[swingTime, stanceTime, stanceLength] = getStepInfo;

high = stanceLength/2;
low = -stanceLength/2;

stepTime = swingDuration+stanceDuration;

timePoints = [0 swingDuration swingDuration+stanceDuration 2*swingDuration+stanceDuration 2*swingDuration+2*stanceDuration];% 3*swingTime+2*stanceTime 3*swingTime+3*stanceTime]; %in second
signalPoints = [low high low high low];% high low];

times = 0:0.0001:(2*swingDuration+2*stanceDuration); %can be finer
times = sort([times timePoints stepTime/6 2*stepTime/6 3*stepTime/6 4*stepTime/6 5*stepTime/6]);
times = unique(times);

oneSignal = interp1(timePoints,signalPoints,times);

I = find(times == stepTime);
I2 = find(times == 2*stepTime);

R1 = oneSignal(1:I-1);%crop to just one cycle.

switch gaitType
    case 'tripod'
        L2 = R1;
        R3 = R1;
        L1Temp = oneSignal(times>stepTime/2);
        L1 = L1Temp(1:length(R1));
        R2 = L1;
        L3 = L1;
    case 'wave1'
        L2Temp = oneSignal(times>=stepTime-stepTime/6);
        L2 = L2Temp(1:length(R1));
        
        R3Temp = oneSignal(times>=stepTime-2*stepTime/6);
        R3 = R3Temp(1:length(R1));
        
        L1Temp = oneSignal(times>=stepTime-3*stepTime/6);
        L1 = L1Temp(1:length(R1));
        
        R2Temp = oneSignal(times>=stepTime-4*stepTime/6);
        R2 = R2Temp(1:length(R1));
        
        L3Temp = oneSignal(times>=stepTime-5*stepTime/6);
        L3 = L3Temp(1:length(R1));
    case 'realTripod' %real tripod is an ideal tripod with realistic delays between within tripod legs.
        L2Temp = oneSignal(times>=stepTime-stepTime*mf_p);
        L2 = L2Temp(1:length(R1));
        
        R3Temp = oneSignal(times>=stepTime-stepTime*hf_p);
        R3 = R3Temp(1:length(R1));
        
        L1Temp = oneSignal(times>=stepTime-stepTime*0.5);
        L1 = L1Temp(1:length(R1));
        
        R2Temp = oneSignal(times>=stepTime-stepTime*(0.5+mf_p));
        R2 = R2Temp(1:length(R1));
        
        L3Temp = oneSignal(times>=stepTime-stepTime*(0.5+hf_p));
        L3 = L3Temp(1:length(R1));
    case 'tetrapod1'
        L2 = R1;
        
        R3Temp = oneSignal(times>=stepTime-stepTime/3);
        R3 = R3Temp(1:length(R1));
        
        L1 = R3;
        
        R2Temp = oneSignal(times>=stepTime-2*stepTime/3);
        R2 = R2Temp(1:length(R1));
        
        L3 = R2;
    case 'tetrapod2'
        L2Temp = oneSignal(times>=stepTime-stepTime/3);
        L2 = L2Temp(1:length(R1));
        
        R3 = L2;
        
        L1Temp = oneSignal(times>=stepTime-2*stepTime/3);
        L1 = L1Temp(1:length(R1));
        
        R2 = L1;
        
        L3 = R1;
end

idealPositions = [R1; R2; R3; L1; L2; L3];

idealTime = times(1:I-1);

idealPostitionsRepeat = repmat(idealPositions,1,11);

idealPhasesRepeat = getPhases(idealPostitionsRepeat);

n=5;
idealPhases = idealPhasesRepeat(:,n*(I-1)+1:(n+1)*(I-1));
end

function phases = getPhases(signals)
phases=zeros(6,size(signals,2));
for i = 1:6
    s=signals(i,:);
    added = 1000;
    oriLength = length(s);
    s = padarray(s, [0 added]);
    
    h = hilbert(s);
    p = angle(h); % Instantaneous phase (wrapped)
    p = p(added+1:oriLength+added);
    phases(i,:) = p;%unwrap(p);
end
end

%% Calculate the score (Compare exp. and ideal phases).
function [relativeExp, patternScore] = calculateGDI(expPhases,idealPhases,expTime,idealTime)
relativeExp=nan;
fakeTime = linspace(expTime(1),expTime(end),length(idealTime));

patternIdealPhases = zeros(6,length(expTime));

%[expPhases, ordered] = unwrapAbout(expPhases,7,[]);
%[idealPhases, ~] = unwrapAbout(idealPhases,7,ordered);
expPhasesW=unwrap(expPhases')';
idealPhasesW=unwrap(idealPhases')';

expPhasesW=expPhasesW-min(min(expPhasesW));
idealPhasesW=idealPhasesW-min(min(idealPhasesW));
for i = 1:6
    patternIdealPhases(i,:) = interp1(fakeTime,idealPhasesW(i,:),expTime);
end

%[pExpT, ~] = unwrapAbout(expPhases,7,[1 3 5]);
%relativeExp = [pExpT(2,:)-pExpT(1,:);pExpT(3,:)-pExpT(1,:);pExpT(4,:)-pExpT(1,:);pExpT(5,:)-pExpT(1,:);pExpT(6,:)-pExpT(1,:)];

% figure
% hold on
% plot(expTime,patternIdealPhases([1 3 5],:),'r');
% plot(expTime,patternIdealPhases([2 4 6],:),'b');
% plot(expTime,expPhases([1 3 5],:),'m');
% plot(expTime,expPhases([2 4 6],:),'c');
% hold off
% close

rms = zeros(1,length(expTime));
iPD=zeros(4,length(expTime));
ePD=zeros(4,length(expTime));
for f = 1:length(expTime)
    
    idealT = patternIdealPhases(:,f);
    expT = expPhasesW(:,f);
  
    idealShift = getDeltas(idealT,'phase','ideal');
    expShift = getDeltas(expT,'phase','exp');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %expShift=wrapNearIdeal(expShift,idealShift,'phase');
    [expShift, idealShift]=wrapExpAndIdeal(expShift,idealShift,'phase');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iPD(:,f)=idealShift;
    ePD(:,f)=expShift;
    
    rms(f) = getGDI(idealShift,expShift);
end
patternScore = rms;

global cycleNum

if false% cycleNum==438 || cycleNum==468
figure;
subplot(5,1,1)
hold on;
plot(expTime,expPhases(1,:),'r');plot(expTime,expPhases(5,:),'g');plot(expTime,expPhases(3,:),'b');plot(expTime,expPhases(4,:),'m');plot(expTime,expPhases(2,:),'k');plot(expTime,expPhases(6,:),'c');
hold off
legend('R1','L2','R3','L1','R2','L3')
ylabel('exp. phase angle')
ylim([-pi pi])

subplot(5,1,2)
hold on;
plot(fakeTime,idealPhases(1,:),'r');plot(fakeTime,idealPhases(5,:),'g');plot(fakeTime,idealPhases(3,:),'b');plot(fakeTime,idealPhases(4,:),'m');plot(fakeTime,idealPhases(2,:),'k');plot(fakeTime,idealPhases(6,:),'c');
hold off
legend('R1','L2','R3','L1','R2','L3')
ylabel('ideal phase angle')
ylim([-pi pi])

subplot(5,1,3)
hold on;
plot(expTime,ePD(1,:),'r');plot(expTime,ePD(2,:),'g');plot(expTime,ePD(3,:),'b');plot(expTime,ePD(4,:),'m')
hold off
legend('L2-R1', 'R3-L2', 'R2-L1', 'L3-R2')
ylabel('exp delta')
ylim([0 pi])

subplot(5,1,4)
hold on;
plot(expTime,iPD(1,:),'r');plot(expTime,iPD(2,:),'g');plot(expTime,iPD(3,:),'b');plot(expTime,iPD(4,:),'m')
hold off
legend('L2-R1', 'R3-L2', 'R2-L1', 'L3-R2')
ylabel('ideal delta')
ylim([0 pi])

subplot(5,1,5)
hold on;
plot(expTime,patternScore,'k')
hold off
xlabel('t')
ylabel('GDI(t)')
ylim([0 7.5])
end
end

function [unwrapedFinal, order] = unwrapAbout(phases,usrCorrect,ordered)
unwrapedFinal=[];
for correctPoint = 1:size(phases,2)
    unwraped = zeros(size(phases,1),size(phases,2));
    for i = 1:size(phases,1)
        unwrapedAfterPoint = unwrap(phases(i,correctPoint:end));
        unwrapedBeforePoint = fliplr(unwrap(fliplr(phases(i,1:correctPoint))));
        unwrapedBeforePoint = unwrapedBeforePoint(1:end-1);
        unwraped(i,:) = [unwrapedBeforePoint unwrapedAfterPoint];
    end
    [~,I] = sort(unwraped(:,correctPoint));
    if isequal(sort(I(4:6))',ordered)
        unwrapedFinal = unwraped;
        %cPnt = correctPoint;
        order = sort(I(4:6))';
        return
    end
    if isempty(ordered)
        if isequal(sort(I(4:6))',[1 3 5]) || isequal(sort(I(4:6))',[2 4 6])
            unwrapedFinal = unwraped;
            %cPnt = correctPoint;
            order = sort(I(4:6))';
            return
        end
    end
end

if isempty(unwrapedFinal)
   unwraped = zeros(size(phases,1),size(phases,2));
    for i = 1:size(phases,1)
        unwrapedAfterPoint = unwrap(phases(i,usrCorrect:end));
        unwrapedBeforePoint = fliplr(unwrap(fliplr(phases(i,1:usrCorrect))));
        unwrapedBeforePoint = unwrapedBeforePoint(1:end-1);
        unwraped(i,:) = [unwrapedBeforePoint unwrapedAfterPoint];
    end
    unwrapedFinal = unwraped;
    order = sort(I(4:6))';
    %cPnt = usrCorrect;
end

end
