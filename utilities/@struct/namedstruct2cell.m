function C = namedstruct2cell(S)
% NAMEDSTRUCT2CELL      Convert struct to cell including field names
% usage: C = namedstruct2cell(S)
%
% Matlab's native struct2cell only gives back the values.  This gives the
% fieldnames and values in a linear cell array: {key1 val1 key2 val2 ...}
%
% 2012-07 phli
%
C = reshape([fieldnames(S) struct2cell(S)]', 1, []);