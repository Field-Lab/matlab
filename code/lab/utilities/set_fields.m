function struc = set_fields(struc, varargin)
% SET_FIELDS    Set multiple struct fields
% usage: struc = set_fields(struc, field_name1, value1, field_name2, value2, ...)
%
% 2010-08 phli
%
for j = 1:2:length(varargin)
    struc.(varargin{j}) = varargin{j+1};
end
