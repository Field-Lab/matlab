function save_monitor_alignment(datarun, the_path)
% save_monitor_alignment     save struct with information needed to align the monitor and array
%
% usage:  datarun = save_monitor_alignment(datarun, the_path)
%
% arguments:     datarun - datarun struct
%               the_path - path to alignment file. if not specified, it will be guessed as
%                           <server_path>/<piece_name>/monitor_alignment.mat
%
%
% 2010-03  gauthier
%



% get information from datarun
var_names = {'points_camera_to_base','T_camera_to_base','T_base_to_array','points_base_to_array','T_base_to_monitor'};
for vv = 1:length(var_names)
    eval(sprintf('%s = datarun.piece.photographic_mapping.%s;',var_names{vv},var_names{vv}))
end

% guess alignment file name if not specified
if ~exist('the_path','var')
    the_path = [server_path guess_piece_name(datarun) '/monitor_alignment.mat'];
end

% save alignment file
save(the_path,var_names{:})

