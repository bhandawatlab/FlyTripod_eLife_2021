%This code allows quick denoising operation on a time series pulse data.
%Similar to low-pass filter.
%@Chanwoo Chun <cc2465@cornell.edu>

function y=denoise(a,n)
%n = mimimum number of consc. 1s that survive
index = [];

fs = find(conv(a,ones(1,n),'valid')==n);

fe = find(conv(fliplr(a),ones(1,n),'valid')==n);

if ~isempty(fs)
    str=fs([true diff(fs)>n]);
    ter = sort(length(a)+1-fe([true diff(fe)>n]));
    index = cell2mat(arrayfun(@(s,e) (s:e), str, ter, 'UniformOutput', false));
end

a = zeros(1,length(a));
a(index) = 1;

y=a;