function loaded = load_enum(enum_class)
% LOAD_ENUM     Load the values (in order) of a Java enum class
% These are 1-indexed for Matlab!
%
% See also: GET_ENUM, LIST_ENUM, FIND_CLASS
%
% 2011-09 phli
%

if ischar(enum_class)
    enum_class = find_class(enum_class);
end

fields = enum_class.getFields();
loaded = cell(length(fields),1);
for i = 1:length(loaded)
    loaded{i} = fields(i).get(enum_class);
end