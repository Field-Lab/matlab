function val = get_gd(h, varargin)
% GET_GD    Get desired (potentially nested) guidata field
% usage: get_gd(handle, [property_names])
%
% 2010-05 phli
%

if nargin < 1
    h = gcf;
end

val = guidata(h);

if ~isempty(varargin)
    val = getfield(val, varargin{:});
end