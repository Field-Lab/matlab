function set_ud(h, varargin)
% SET_UD    Set a field of the UserData property, assuming UD is a struct
% usage: set_ud(handle, property_name[s], value)
% 
% NOTE: UserData can probably be retired in favor of SETAPPDATA/GETAPPDATA
%
% In general, SNL-E policy should be to keep UD setup as a struct.  SET_UD
% confirms that this is the case and then sets the desired (potentially
% nested) field.
%
% phli 2010-05
%


ud = get(h, 'UserData');

if isempty(ud)
    ud = struct();
end

if ~isstruct(ud)
    warning('SET_UD:UserDataNotStruct', 'UserData is already set and not a struct for this handle');
    return;
end

ud = setfield(ud, varargin{:});
set(h, 'UserData', ud);