function deltas = getDeltas(A,flag1,flag2)
%R1[1] R2[2] R3[3] L1[4] L2[5] L3[6]
%
%L2[5]-L1[4], L3[6]-L2[5], R2[2]-R1[1], R3[3]-R2[2],
%R3[3]-L1[4]+L3[6]-R1[1]
%d1 = A(5)-A(4); d2 = A(6)-A(5); d3 = A(2)-A(1); d4 = A(3)-A(2);
%d5 = A(3)-A(4)+A(6)-A(1);
%deltas = [d1 d2 d3 d4 d5];

%L2[5]-R1[1], R3[3]-L2[5], R2[2]-L1[4], L3[6]-R2[2]
d1 = A(5)-A(1); d2 = A(3)-A(5); d3 = A(2)-A(4); d4 = A(6)-A(2);
deltas = [d1 d2 d3 d4];

if strcmp(flag1,'phase')
   deltas = -deltas; 
elseif strcmp(flag1,'time')
    %no need to do anyting.
end

end