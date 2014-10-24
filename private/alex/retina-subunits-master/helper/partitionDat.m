function [train test] = partitionDat(dat,trainFrac,nBins,dim)

% [train test] = partitionDat(dat,trainFrac,nBins,dim)
%
% subdivision a data set into chunks (# = nBins), from each chunk, grab a
% fraction for training and a fraction for testing (testFrac);
% dim is the dimension along which to partition

if nargin < 3
    nBins = 10;
    dim = 1;
end

if nargin < 4
    dim = 1;
end

len = size(dat,1);
dBin = len/nBins;
binEdges = round([1:dBin:len]);
binEdges(end+1) = len + 1;

train = [];
test = [];
for i=1:nBins
    tmpWidth = binEdges(i+1)-binEdges(i);
    if dim == 1
        train = [train; dat(binEdges(i):binEdges(i)+round(trainFrac*tmpWidth),:)];
        test = [test; dat(binEdges(i)+round(trainFrac*tmpWidth)+1:binEdges(i+1)-1,:)];
    end
end



