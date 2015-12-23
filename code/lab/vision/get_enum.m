function enum = get_enum(enum_class, indexorname)
% GET_ENUM      Get the Java object corresponding to the given value from the given ENUM_CLASS.
% The value can be selected as a numerical index (1-indexed for Matlab) 
% into the enum, or else as a string of the name of the value.
%
% See also: FIND_CLASS, LIST_ENUM, LOAD_ENUM
%
% 2011-09 phli
%

if ischar(enum_class)
    enum_class = find_class(enum_class);
end

if isnumeric(indexorname)
    fields = enum_class.getFields();
    enum = fields(indexorname).get(enum_class);
elseif ischar(indexorname)
    enum = java.lang.Enum.valueOf(enum_class, indexorname);
else
    enum = [];
end