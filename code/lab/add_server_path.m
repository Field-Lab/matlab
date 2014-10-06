function full_path = add_server_path(input_path)
% ADD_SERVER_PATH    Add server path, checking to make sure path is not already absolute
% usage: full_path = add_server_path(input_path)
%
% 2010-05 phli, abstracted out of load_data/expand_rrs_path, added windows compatibility
% 2010-08 phli, now simply calls abstracted version ADD_BASE_PATH; most calls 
%               to this should probably be replaced with add_base_path(inputpath, server_path());
%

full_path = add_base_path(input_path, server_path());