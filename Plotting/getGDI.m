function gdi = getGDI(u,v)
%u=u/norm(u);
%v=v/norm(v);
%gdi = sum(u.*v);
gdi=sqrt(sum((u-v).^2));
end