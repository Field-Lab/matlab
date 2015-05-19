function str_cell = join(str, delim)
% Join a cell array of strings into a single char array, separated by a given deliminter
%
% 2010-01 phli

% ToDo This can probably be made more efficient
delims = repmat({delim}, size(str(:)'));
joined = vertcat(str(:)', delims);
str_cell = [joined{1:end-1}];