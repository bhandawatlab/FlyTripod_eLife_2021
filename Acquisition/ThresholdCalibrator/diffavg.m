%Order of operation in this function:
% 1) take 3D matrix
% 2) take absolute "diff" of the matrix in third dimension
% 3) sum all the elements in each thrid dimension (output=1D)
% 4) take average value of elements in the 1D matrix
% 5) multiply a predefined factor
%This quantifies how much "motion" there was in a given video.
%
% @Chanwoo Chun, <cc2465@cornell>

function y = diffavg(F)

n = size(F,3);

diff = zeros(size(F,1), size(F,2), n-1);
tdiff = zeros(n-1,1);

for fn = 1:n-1
    diff(:,:,fn) = abs(F(:,:,fn)-F(:,:,fn+1));
    tdiff(fn)= sum(sum(diff(:,:,fn),1),2);
end

tdiffavg = mean(tdiff)*0.00001;

y=tdiffavg;



