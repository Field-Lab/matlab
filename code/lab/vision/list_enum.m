function list = list_enum(enum_class)
% LIST_ENUM     Show the possible values (in order) of a Java enum class
% These are 1-indexed for Matlab!
%
% See also: LOAD_ENUM, GET_ENUM, FIND_CLASS
%
% 2011-09 phli
%

if ischar(enum_class)
    enum_class = find_class(enum_class);
end

fields = enum_class.getFields();
list = cell(length(fields),1);
for i = 1:length(list)
    list{i} = fields(i).getName().toCharArray()';
end