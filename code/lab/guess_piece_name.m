function piece_name = guess_piece_name(datarun)
% guess_piece_name     Identify the name of the piece
%
% usage:  piece_name = guess_piece_name(datarun)
%
% arguments:     datarun - datarun struct
%
% outputs:     piece_name - string, guessed name
%
%
% e.g. if datarun has field datarun.names.rrs_prefix = '/snle/lab/Experiments/Array/Analysis/2005-07-06-2/data000-gdf/data000'
%       then guess_piece_name(datarun) -->  '2005-07-06-2'
%
%
%
% 2010-01  gauthier
%




% GET PREFIX

if isfield(datarun,'names')
    % get rrs_prefix, if it exists
    if isfield(datarun.names,'rrs_prefix')
        prefix = datarun.names.rrs_prefix;
    elseif isfield(datarun.names,'rrs_neurons_path')
        % otherwise, use neurons file prefix
        prefix = datarun.names.rrs_neurons_path;
    end
end

if ~exist('prefix','var')
    piece_name = '';
    return
end





% CHECK PREFIX FOR DATE


% check whether the server path is in the prefix
[junk,server_path_end] = regexp(prefix,server_path);

% if so use it to find piece name
if ~isempty(server_path_end)
    
    % get where piece name starts
    piece_name_start = server_path_end + 1;
    
    % assume it ends at the next '/'
    temp = regexp(prefix(piece_name_start:end),'/');
    
    % if there is a '/'...
    if ~isempty(temp)
        piece_name_end = piece_name_start+temp(1)-2;
        
        % extract the piece name
        piece_name = prefix(piece_name_start:piece_name_end);
        return
    else
        % otherwise, I have no idea what's going on
        piece_name = '';
        return
    end
    
else
    % if file is not on the server...
    fprintf('\nCurrently ''guess_piece_name'' only works for dataruns on the server.  Perhaps YOU can expand its functionality!\n\n')
    piece_name = '';
    return
    
    % the code here should use regular expression to check for a date, like "20##-##-##-#" or "19##-##-##-#"
end


    
    