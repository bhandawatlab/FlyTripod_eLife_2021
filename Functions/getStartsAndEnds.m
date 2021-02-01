function [starts, ends] = getStartsAndEnds(v)
    if iscolumn(v)
        v = v';
    end
    w = [false v~=0 false]; %// "close" v with zeros, and transform to logical
    starts = find(w(2:end) & ~w(1:end-1)); %// find starts of runs of non-zeros
    ends = find(~w(2:end) & w(1:end-1))-1; %// find ends of runs of non-zeros
end