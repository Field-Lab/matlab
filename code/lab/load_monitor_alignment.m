function datarun = load_monitor_alignment(datarun, the_path)
% load_monitor_alignment     load struct storing information needed to align the monitor and array
%
% usage:  datarun = load_monitor_alignment(datarun, the_path)
%
% arguments:     datarun - datarun struct with field
%               the_path - path to alignment file. if not specified, it will be guessed as
%                           <server_path>/<piece_name>/monitor_alignment.mat
%
% outputs:     datarun - datarun struct with fields as explained in compute_monitor_to_array_transformation
%
%
% 2010-03  gauthier
%




% guess alignment file name if not specified
if ~exist('the_path','var')
    the_path = [server_path guess_piece_name(datarun) '/monitor_alignment.mat'];
end

% load alignment file
S = load(the_path);

% put the information in datarun
if ~isfield(datarun.piece, 'photographic_mapping')
    datarun.piece.photographic_mapping = struct();
end
datarun.piece.photographic_mapping = structmerge(datarun.piece.photographic_mapping, S);

% generate the rest of it
datarun = compute_monitor_to_array_transformation(datarun, 'reselect', false);



function structout = structmerge(structin, newfields)
% Set the output to the input
structout = structin;

% Loop over all the fields of newfields and assign them to the output
fields = fieldnames(newfields);
for i = 1:length(fields),
    structout.(fields{i}) = newfields.(fields{i});
end