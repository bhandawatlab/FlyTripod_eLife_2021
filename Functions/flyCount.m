function flyNumCount=flyCount(flyNums)

flyNumCount=NaN(max(flyNums),1);
for i=1:length(flyNumCount)
    flyNumCount(i)=length(find(flyNums==i));
end

end