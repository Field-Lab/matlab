function str_cell = split(str, delim)
% SPLIT    Break a char array into substrings at given delimiter
% example: split('foo:bar', ':')
%          >> {'foo', 'bar'} 
%
% 2010-01 phli
%

delimlen = length(delim);
str_cell = {};

i = strfind(str, delim);
while ~isempty(i)
    str_cell{end+1} = str(1:i(1)-1); %#ok<AGROW>
    str = str(i(1)+delimlen:end);
    i = strfind(str, delim);
end
str_cell{end+1} = str;