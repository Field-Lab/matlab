function datarun = expand_rrs_prefix(datarun, varargin)
% EXPAND_RRS_PREFIX     expand rrs_prefix to the names for any rrs files that exist
%
% usage:  datarun = expand_rrs_prefix(datarun)
%
% arguments:  datarun - datarun struct with field datarun.names.rrs_prefix
%
% outputs:    datarun - datarun struct with these fields filled in if files exist
%
%       datarun.names.rrs_neurons_path
%       datarun.names.rrs_params_path
%       datarun.names.rrs_ei_path
%       datarun.names.rrs_sta_path
%
%   NOTE: existing values are NOT over-written
%
%
%  A message is printed if no new fields were filled in from the prefix.
%
% gauthier 2008-03
% phli     2011-01, Added cov files
% phli     2011-10, Added split files
% phli     2012-07, Added verbose option via inputParser
%

% ensure rrs_prefix exists
if ~isfield(datarun,'names') || ~isfield(datarun.names,'rrs_prefix'), return; end

opts = inputParser();
opts.addParamValue('quiet', false);
opts.parse(varargin{:});
opts = opts.Results;

% fields to check
rrs_fields = {'rrs_neurons_path','rrs_params_path','rrs_ei_path','rrs_sta_path','rrs_globals_path','rrs_movie_path','rrs_cov_path','rrs_ncov_path','rrs_wcov_path'};

% if the prefix does not start with a '/', add server path to the beginning
datarun.names.rrs_prefix = add_server_path(datarun.names.rrs_prefix);


% initialize
no_files_found = true;

% go through each field name
for rr = 1:length(rrs_fields)
    
    % don't overwrite existing field, unless it's empty
    if ~isfield(datarun.names,rrs_fields{rr}) || isempty(datarun.names.(rrs_fields{rr}))
        
        % generate the name of the file
        prospective_file = [datarun.names.rrs_prefix '.' strrep(strrep(rrs_fields{rr},'_path',''),'rrs_','')];
        
        % if it exists, save it in datarun.names
        if exist(prospective_file,'file')
            datarun.names.(rrs_fields{rr}) = prospective_file;
            
            % note that a file was found
            no_files_found = false;
        end
        
        % split files
        if exist([prospective_file '-0'], 'file')
            datarun.names.(['split_' rrs_fields{rr}]) = prospective_file;
            
            % note that a file was found
            no_files_found = false;
        end        
    end
end

% if no files were found, print message
if no_files_found && ~opts.quiet
    fprintf('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nNOTE: No files added from this prefix: %s.*\n   Either the path is incorrect, or all files were already specified in datarun.names\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n',...
        datarun.names.rrs_prefix)
end