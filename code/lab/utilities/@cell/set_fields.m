function ca = set_fields(ca, varargin)
% SET_FIELDS    Set fields for each element of a cell array of structs
% usage: cell_array = set_fields(cell_array, field_name1, value1, field_name2, value2, ...)
%
% 2010-08 phli
%

for i = 1:numel(ca)
    elem = ca{i};

    if isstruct(elem)
        elem = set_fields(elem, varargin{:});
        ca{i} = elem;
    end
end