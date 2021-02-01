function [wrapped_expD, wrapped_idealD] = wrapExpAndIdeal(expD,idealD,flag)

if strcmp(flag,'time')
    wrapped_expD=wrapToPi(expD*2*pi)/(2*pi);
    wrapped_idealD=wrapToPi(idealD*2*pi)/(2*pi);
else
    wrapped_expD=wrapToPi(expD);
    wrapped_idealD=wrapToPi(idealD);
end

end