%This code returns a reasonably good background estimation of a video.
%@Chanwoo Chun <cc2465@cornell.edu>

function background = trueBackground(frames)
tframes=size(frames,3);

for column = 1:size(frames,2)
    %A single column of a frame is chosen. Now, from the start to the end
    %of the video, take a sum of the pixel intensity values of the column
    %and then square it.
    for i = 1:tframes
        columnBrSum(i) = sum(frames(:,column,i))^2;
    end
    %Imagin you are plotting columnBrSum array (x-index, y-element
    %value).
    %If some object passed through this column at some point, the array
    %columnBrSum will have some elements whose values are higher than
    %others.
    %Therefore, when columnBrSum is plotted, there would be a baseline
    %where nothing happens. We want to identify the index values of the
    %baseline, and take average of the columnes of that index (frame)
    %values. Doing this for every column will give us a background.
    
    
    %Again, imagine the columnBrSum plot. Now, we add a horizontal line
    %that crosses on top of the columnBrSum plot. This line starts from the
    %minimum value of columnBrSum and slides up a little by a predefined
    %amount (this value may need to be changed depending on a
    %characteristic of a video). In the for loop below, we are essentially
    %counting the number of intersections between the line and columnBrSum
    %plot. 
    lowMinimum = min(columnBrSum);
    maxPopNum = 1;
    for scan = lowMinimum:0.05e7:lowMinimum+5e7
        population = false(tframes,1);
        for i = 2:tframes
            if (columnBrSum(i)>=scan && columnBrSum(i-1)<=scan)||...
                    (columnBrSum(i)<=scan && columnBrSum(i-1)>=scan)
                population(i-1:i) = true;
            end
        end

        popNum = size(find(population),1);

        if popNum > maxPopNum
            maxPop = population;
            maxPopNum = popNum;
        end
    end

    
    %Now the following block determines when the object appeared and
    %disappeared. We are going to mask out the portion where the object is
    %detected, and then take average of the columns.
    %However, other simpler methods may work. For example, we can skip to
    %the line where variable multiBack is defined, and then set multiBack
    %equal to population <multiBack=population;>. For this, the below 
    %for-loop can be commented out.
    falseCount = 0;
    maxCount = 0;
    for i = 1:tframes   
        if maxPop(i)==false
            falseCount = falseCount+1;
            if i~=tframes
                if maxPop(i+1)==true
                   if  falseCount > maxCount
                       maxCount = falseCount;
                       gone = i;
                       appear = gone-maxCount;
                   end
                   falseCount = 0;
                end
            else
                if  falseCount > maxCount
                       maxCount = falseCount;
                       gone = i;
                       appear = gone-maxCount;
                end
            end
        end
    end
    
    if appear == 0
       appear = 1;
    end

    mask = true(tframes,1);
    mask(appear:gone) = false;

    multiBack = squeeze(frames(:,column,mask));

    avgBackC(:,column) = mean(multiBack,2);
    
    appear = 0;
    gone = 0;
end

background = uint8(255*mat2gray(avgBackC));
