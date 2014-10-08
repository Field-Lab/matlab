function extractedData = extractDataSubset3D(data, d1Extract, d2Extract, d3Extract)

% This function extracts a subarray from within the 3D array "data" defined by d1Extract, d2Extract,
% and d3Extract (denoting dimensions 1, 2 and 3).  Only elements corresponding to a value of 1 in
% each d#Extract are extracted, and the rest of the elements are left out.  
% Extracted elements stay in the same order.
%
% arguments
%   data: a 3D array
%   d#Extract: 1D binary arrays denoting which data to extract (extract = 1, don't extract = 0)
%
% returns
%   extractedData: the subset of elements data(i,j,k) for which d1Extract(i) = d2Extract(j) =
%   d3Extract(k) = 1


if length(d1Extract) ~= size(data, 1)
    error('The length of "d1Extract" is not equal to the size of the first dimension of "data."  Aborting.')
end

if length(d2Extract) ~= size(data, 2)
    error('The length of "d2Extract" is not equal to the size of the second dimension of "data."  Aborting.')
end

if length(d3Extract) ~= size(data, 3)
    error('The length of "d3Extract" is not equal to the size of the third dimension of "data."  Aborting.')
end

for i = -size(data, 1):-1
    if ~d1Extract(-i)
        data(-i,:,:) = [];
    end
end

for i = -size(data, 2):-1
    if ~d2Extract(-i)
        data(:,-i,:) = [];
    end
end

for i = -size(data, 3):-1
    if ~d3Extract(-i)
        data(:,:,-i) = [];
    end
end

extractedData = data;