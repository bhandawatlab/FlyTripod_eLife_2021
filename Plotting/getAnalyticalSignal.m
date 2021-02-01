%This function extracts analytical signal of leg position through hilbert
%transform.
%This function is called from cycle_analysis.m
%
%@Chanwoo Chun, <cc2465@cornell.edu> Jan. 18, 2019

function analyticalSignals=getAnalyticalSignal(legPosStruct,vStarts,vEnds)
legPosY(1,:)=legPosStruct.R1(:,2)';
legPosY(2,:)=legPosStruct.R2(:,2)';
legPosY(3,:)=legPosStruct.R3(:,2)';
legPosY(4,:)=legPosStruct.L1(:,2)';
legPosY(5,:)=legPosStruct.L2(:,2)';
legPosY(6,:)=legPosStruct.L3(:,2)';

allPhases = NaN(size(legPosY,1),size(legPosY,2));
%Get analytical signal in each valid section.
for i = 1:length(vStarts)
   %extract one valid section
   vLegPos=legPosY(:,vStarts(i):vEnds(i));
   phases=NaN(6,size(vLegPos,2));
   for legNum = 1:6
      %extract one leg
      stancePos=vLegPos(legNum,:);
      %find stance start and end times.
      %We first need to interpolate the swing phases.
      [starts, ends] = getStartsAndEnds(stancePos);
      if isempty(ends)
          continue
      end
      %loop through stance ends.
      for j = 1:length(ends)
         currentStanceEnd = ends(j); 
         nextStanceStart = min(starts(starts>currentStanceEnd));
         %for the last stance end, there could be no stance start after
         %that. In that case, omit the iteration.
         if isempty(nextStanceStart)
             continue
         end
         currentPos = stancePos(currentStanceEnd);
         nexPos = stancePos(nextStanceStart);
         numData = nextStanceStart-currentStanceEnd+1;
         %interpolate swing phase
         swingPos = interp1([1 numData],[currentPos nexPos],1:numData);
         
         %remove last and first frame because those are for stance phases.
         swingPos(end)=[];
         swingPos(1)=[];
         stancePos(currentStanceEnd+1:nextStanceStart-1)=swingPos;
      end
      cumsummed=cumsum(stancePos).*cumsum(stancePos,'reverse');
      stancePos(cumsummed==0)=nan;
      stancePos = stancePos - (max(stancePos,[],'omitnan')+min(stancePos,[],'omitnan'))/2;
      stancePos(cumsummed==0)=0;
      
      %Remember, this phase may include swing phases that couldn't be
      %interpolated, since they are at the beginning or end of the time
      %series. They were left as zeros. So do not believe the phase values
      %at the beginning or end. To fix this, we need the next step.
      phase = getPhases(stancePos);
      

      %Now we have leg position for stance and swing. However, if the time
      %series started with swing phase, the array will have one or more
      %zero(s) at the beginning. Also, if the time series ended with swing
      %phase, the array will have one or more zero(s) at the end. We will
      %replace these zeros will NaN. Once we do this for all the legs, we
      %will omit any time where at least one leg is in NaN state. We will
      %do that because we cannot get relative phase difference between the
      %legs, if one of the legs being compared is in NaN state.
      phase(cumsummed==0)=NaN;
      phases(legNum,:)=phase;
   end
   allPhases(:,vStarts(i):vEnds(i))=phases;
end

for legNum = 1:6
allPhases(legNum,isnan(sum(allPhases,1)))=NaN;
end

analyticalSignals=allPhases;
end

%% Get Phases of all six legs.
function phase = getPhases(signal)
    added = 1000;
    oriLength = length(signal);
    signal = padarray(signal, [0 added]);
    
    h = hilbert(signal);
    p = angle(h); % Instantaneous phase (wrapped)
    phase = p(added+1:oriLength+added);

end