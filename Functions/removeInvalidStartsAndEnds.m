function [cleaned_starts, cleaned_ends]=removeInvalidStartsAndEnds(starts,ends,validStarts,validEnds,option)
        %option1: Remove starting and ending indices if they conincide with
        %starting and ending indices of validity.
        %option2: Crop out the first and last step if they are incomplete.
        if strcmp(option,'option1')
        [cleaned_starts, ~] = setdiff(starts,validStarts); %remove starting indices that are the same as valid start indices
        [cleaned_ends,~] = setdiff(ends,validEnds); %same for the ending indices
        elseif strcmp(option,'option2')
        [starts, iS] = setdiff(starts,validStarts); %remove starting indices that are the same as valid start indices
        ends = ends(iS); %remove corresponding ending indices FOR DURATION CALCULATION
        [cleaned_ends,iE] = setdiff(ends,validEnds); %same for the ending indices
        cleaned_starts = starts(iE); %same FOR DURATION CALCULATION
        end
end