function val = get_ud(h, varargin)
% GET_UD    Get desired UserData field, assuming UD is a struct
% usage: get_ud(handle, property_name[s])
%
% In general, SNL-E policy should be to keep UserData setup as a struct.
% This function will retrieve the desired (potentially nested) field from
% UD.
%
% 2010-05 phli
%

if nargin < 1
    h = gcf;
end

val = get(h, 'UserData');

if ~isempty(varargin)
    val = getfield(val, varargin{:});
end