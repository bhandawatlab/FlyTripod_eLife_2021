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
