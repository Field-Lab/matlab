function full_path = add_base_path(input_path, base_path)
% ADD_BASE_PATH    Add server path, checking to make sure path is not already absolute
% usage: full_path = add_base_path(input_path, base_path)
%
% 2010-05 phli, abstracted out of load_server_path
%

if ~strcmp(input_path(1),'/') && ~is_windows_root_path(input_path)
    full_path = fullfile(base_path, input_path);
else
    full_path = fullfile('',input_path);
end