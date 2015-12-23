function datarun = load_index(datarun,varargin)
% LOAD_INDEX     Load information about a series from an index file
%
% usage:  datarun = load_index(datarun)
%
% arguments:  datarun - datarun struct with fields specifying the index file and series
%                         (datarun.names.experiment - piece number, e.g. '2005-04-26-1'
%                           datarun.names.condition - experiment condition, e.g. 'rf-6-bk1-jg-2')
%              params - struct of optional parameters (see below)
%
% outputs:    datarun - datarun struct with the following fields added, as possible
%
%   file information (existing values are NOT over-written, unless they are empty)
%
%       datarun.names.rrs_prefix
%       datarun.names.rrs_neurons_path
%       datarun.names.rrs_params_path
%       datarun.names.rrs_ei_path
%       datarun.names.rrs_sta_path
%       datarun.names.movie_xml_path
%       datarun.names.nickname
%
%
% 	piece information (existing values are over-written)
%
%       datarun.piece.eye
%       datarun.piece.eccentricity_radius
%       datarun.piece.eccentricity_clock
%       datarun.piece.rig
%      	datarun.piece.preparation
%       datarun.piece.optical_path_direction
%       datarun.piece.objective
%       datarun.piece.recording_system
%       datarun.piece.alive_montage_lores
%       datarun.piece.alive_montage_hires
%       datarun.piece.photographic_mapping_pm2
%       datarun.piece.photographic_mapping_pm10
%       datarun.piece.photographic_mapping_pm32
%       datarun.piece.photographic_mapping_array
%       datarun.piece.photographic_mapping_array_edges
%
%       datarun.piece.array *** not implemented yet
%     	datarun.piece.aperture *** not implemented yet
%
%
%   cell type information (existing values are over-written)
%
%       datarun.obvius.cell_types
%
%
%   stacks initialized in datarun.stacks
%       array coordinates placeholder is initialized automatically if stacks doesn't already exist
%       array image is initialized if datarun.piece.array_id is found
%       alive_montage and photographic_mapping fields are used to populate appropriate stacks
%
%
%
% optional params, their default values, and what they specify:
%
% load_piece_info        	true        	load piece information
% load_rrs_paths        	true        	load rrs paths
% load_and_sync_stimulus  	true        	load and stimulus, and sync with datarun.stimulus
% load_cell_types        	true        	load cell types
% index_path                []              allows entering an alternative index path
%
%
%
%
% 2006-03-16    greschner, shlens
% 2006-07-24    shlens
% 2008-03       gauthier
% 2009-11       gauthier, added code to load piece info and stimulus
% 2010-01       phli, added code to load cell_file from Index, added ')' stripping to convert_value
% 2010-03       phli, added check for pseudo-condition based on rrs_prefix, e.g. looks for :condition :data001.  Added ')' stripping to convert_path
% 2010-03       gauthier, added ability to correctly read quoted items with spaces in them, e.g. file names
% 2010-09       phli, migrating to new stacks system for storing images
% 2014-05       gdf, added the ability to specify arbitary path to index file
%

% TO BE IMPLEMENTED
%
% * read stimulus information
%
%




% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('load_piece_info',        true);
p.addParamValue('load_rrs_paths',         true);
p.addParamValue('load_and_sync_stimulus', true);
p.addParamValue('load_cell_types',        true);
p.addParamValue('cell_type_overwrite',    false);
p.addParamValue('basepath',               server_data_path);
p.addParamValue('index_path', []);


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




%  READ IN ENTIRE INDEX FILE

% get piece name
if isfield(datarun,'names') && isfield(datarun.names,'experiment')
    piece_name = datarun.names.experiment;
else
    % guess if needed
    piece_name = guess_piece_name(datarun);
end

% index file location
if isempty(p.Results.index_path)
    index_path = [server_path piece_name '/Index'];
else
    index_path = p.Results.index_path;
end
    
% open index file
fid=fopen(index_path);
if (fid == -1)
    error('load_index: can''t find file %s',index_path);
end

% read entire Index file
file=textscan(fid, '%s','Whitespace',' \b\t\n');
fclose(fid);
index_text=file{1};

% remove white space
for ii=1:length(index_text)
    index_text{ii} = strtrim(index_text{ii});
end




% LOAD PIECE INFORMATION
if params.load_piece_info

    % find text related to the piece
    piece_text = get_piece_text(index_text,index_path);

    % make list of fields to save
    piece_fields = {'eye','eccentricity_radius','eccentricity_clock','rig','recording_system',...
        'preparation','optical_path_direction','objective','array_id','array',...
        'alive_montage_lores','alive_montage_hires'...
        'photographic_mapping_pm2','photographic_mapping_pm10','photographic_mapping_pm32',...
        'photographic_mapping_array','photographic_mapping_array_edges', 'display'};

    % find these fields in the index file
    % note: index_paths.X = [] if not specified in index file
    piece_info = get_rrs_fields(piece_text, piece_fields);
    
    
    % handle array map specially
    if ~isempty(piece_info.array)
        if mod(length(piece_info.array),2)==0
            piece_info.array = reshape(convert_value(piece_info.array)',2,[])';
        else
            error('array field must have an even number of elements, but it has %d',length(piece_info.array))
        end
    end
        

    % go through each field name
    for rr = 1:length(piece_fields)
        if isnumeric(piece_info.(piece_fields{rr}))
            datarun.piece.(piece_fields{rr}) = piece_info.(piece_fields{rr});
        else
            % and save it in datarun.piece
            datarun.piece.(piece_fields{rr}) = convert_value(piece_info.(piece_fields{rr}));
        end
    end

    
    % Setup images to work with new stacks system
    datarun = copy_images_to_stacks(datarun, 'basepath', params.basepath);
end




% Look for condition info, if exists, otherwise end here
if isfield(datarun,'names') && isfield(datarun.names,'condition')

    % Get the text describing this condition
    condition_text = get_condition_text(index_text,datarun.names.condition);
    
    % give an error if the condition wasn't found
    if isempty(condition_text)
        error('load_index: condition ''%s'' not found in index file %s .',datarun.names.condition,index_path)
    end

elseif isfield(datarun, 'names') && isfield(datarun.names, 'rrs_prefix')

    % Try looking for a pseudo-condition, e.g. 'data001'
    rrs_prefix_split = split(datarun.names.rrs_prefix, '/');
    dataname = rrs_prefix_split{end};
    condition_text = get_condition_text(index_text, dataname);
    
    if isempty(condition_text)
        return
    end
    
else
    return
end






% LOAD INFORMATION FOR datarun.names
% don't overwrite existing values
if params.load_rrs_paths

    % set list of rrs fields
    rrs_fields = {'rrs_prefix','rrs_neurons_path','rrs_params_path','rrs_ei_path','rrs_sta_path','nickname','movie_xml','movie_xml_path'};

    % get rrs fields from index file
    % note: index_paths.X = [] if not specified in index file
    index_paths = get_rrs_fields(condition_text, rrs_fields);

    % go through each field name
    for rr = 1:length(rrs_fields)
        % don't overwrite existing field, unless it's empty
        if ~isfield(datarun.names,rrs_fields{rr}) || isempty(datarun.names.(rrs_fields{rr}))
            % load path from index file if it was specified
            if ~isempty(index_paths.(rrs_fields{rr}))
                % for some fields, just copy the value over
                if strcmpi(rrs_fields{rr},'nickname')
                    datarun.names.(rrs_fields{rr}) = index_paths.(rrs_fields{rr});
                else % for most fields, treat the value as a path that must be converted
                    datarun.names.(rrs_fields{rr}) = convert_path(index_paths.(rrs_fields{rr}), rrs_fields{rr});
                end
            end
        end
    end
end


% Get cell_file, if not already filled
if ~isfield(datarun.names, 'cell_file') || isempty(datarun.names.cell_file)
    index_paths = get_rrs_fields(condition_text, {'cell_file'});
    if ~isempty(index_paths.cell_file)
        datarun.names.cell_file = index_paths.cell_file;
    end
end


% LOAD INFORMATION FOR datarun.stimulus
% sync stimulus values here with the existing ones in datarun
if params.load_and_sync_stimulus

    % set list of rrs fields
    stimulus_fields = {'independent','stixel_height','stixel_width',...
        'field_height','field_width','interval','x_start','x_end','y_start','y_end',...
        'monitor_x','monitor_y'};

    % get rrs fields from index file
    % note: index_paths.X = [] if not specified in index file
    stimulus = get_rrs_fields(condition_text, stimulus_fields);

    % translate values to matlab format, and remove empty fields
    for rr = 1:length(stimulus_fields)
        if ~isempty( stimulus.(stimulus_fields{rr}) )
            stimulus.(stimulus_fields{rr}) = convert_value(stimulus.(stimulus_fields{rr}));
        else
            stimulus = rmfield(stimulus,stimulus_fields{rr});
        end
    end

    % sync with existing stimulus
    datarun = sync_stimulus(datarun,stimulus,sprintf('index file ''%s''',index_path));

end




% LOAD CELL TYPES
if params.load_cell_types

    % from the previous code on the server

    % load cell types
    %cell_types = load_cell_types(index_path, datarun.names.condition);
    cell_types = get_cell_types(condition_text);

    % if any were found...
    if ~isempty(fieldnames(cell_types))

        % arrange in new format
        f=fieldnames(cell_types);
        for i=1:length(f)
            cell_types_cell{i}.name = f{i};
            cell_types_cell{i}.cell_ids = getfield(cell_types,f{i})';
        end

        % order cell types, and save them in datarun
        datarun.obvius.cell_types = order_cell_types(cell_types_cell);

        % load these cell types is none have been defined so far
        if ~isfield(datarun,'cell_types') || isempty(datarun.cell_types) || cell_type_overwrite
            datarun.cell_types = datarun.obvius.cell_types;
        end
    else
        % if no cell types were found, make an empty struct, but don't overwrite previous
        if ~isfield(datarun,'obvius') || ~isfield(datarun.obvius,'cell_types') || isempty(datarun.obvius.cell_types)
            datarun.obvius.cell_types = struct;
        end
    end


end







% HELPER FUNCTIONS


function piece_text = get_piece_text(index_text,index_path)
% get text from the piece info

% scan until finding first close parenthesis
for ii = 1:length(index_text)
    if strcmpi(index_text{ii},')')
        % end piece info text here
        piece_info_end = ii;
        break
    end
end

% take text between them
if exist('piece_info_end','var')
    piece_text = {index_text{1:piece_info_end-1}};
else
    error(sprintf('Can not find end of piece info for file %s\nBe sure the closing parenthesis is preceded by white space.',index_path))
end



function condition_text = get_condition_text(index_text,condition)
% get text from the index file which corresponds to the desired condition

% initialize variable to store the condition start location
% once nonempty, this indicates that the start location has been found
condition_start = [];

% scan until the right condition is found
for ii = 1:length(index_text)

    % if word ii matches the desired condition
    if (strcmpi(index_text{ii},condition));
        % check to be sure the previous word was 'condition'
        if regexpi(index_text{ii - 1},'condition')
            % if so, note that this is where the condition begins
            condition_start = ii - 1;
            continue
        end
    end

    % if the beginning of the desired condition has been found
    if ~isempty(condition_start)
        % search for its end (by looking for the beginning of the next condition)
        if regexpi(index_text{ii},'condition')
            condition_end = ii - 1;
            break
        end
    end
end

% if the condition could not be found, return an empty array
if isempty(condition_start)
    condition_text = [];

else % otherwise, return the words between the start and end points
    if exist('condition_end','var')
        % if another condition followed, end there
        condition_text = {index_text{condition_start:condition_end}};
    else
        % if there was no condition after the desired condition, assume it was the last condition
        % in that case, use the end of the index file as the endpoint
        condition_text = {index_text{condition_start:end}};
    end
end



function results = get_rrs_fields(condition_text, desired_fields)
% returns field values from the condition text with some character substitution


% in the index file, dashes are used in place of underscores, so search for these names
desired_field_names = strrep(desired_fields,'_','-');

% go through each field name
for dd = 1:length(desired_fields)
    % scan the condition text until finding the field name
    for ii = 1:length(condition_text)
        if (strcmpi(condition_text{ii},[':' desired_field_names{dd}]));
            
            
            
            % HANDLE QUOTES
            % if the value is in quotes and spans multiple items, be sure to grab the whole thing
            
            if condition_text{ii+1}(1) == '"' && condition_text{ii+1}(end) ~= '"'
                
                % find the closing quotes
                jj = ii+1;
                % cycle through all remaining items
                while jj < length(condition_text)
                    jj = jj +1;
                    % if it begins with quotes, give error
                    if condition_text{jj}(1) == '"'
                        error('could not resolve quoted item starting with << %s >> because a subsequent item started with "',condition_text{ii+1})
                        % if it ends with quotes, proceed
                    elseif condition_text{jj}(end) == '"'
                        break
                    end
                    % otherwise keep going
                end
                % if you go to the end without finding quotes, give error
                if condition_text{jj}(end) ~= '"'
                    error('could not resolve quoted item starting with << %s >> because no subsequent " were found',condition_text{ii+1})
                end
                % otherwise, concatenate the words with spaces inserted
                field_value = [];
                for kk=ii+1:jj
                    field_value = [field_value ' ' condition_text{kk}];
                end
                % strip the space at the begining
                field_value = field_value(2:end);
                
                
                
                
            % HANDLE PARENTHESIS
            % if the value starts with open parenthesis, be sure to grab the whole thing

            elseif condition_text{ii+1}(1) == '('
                
                % initialize variable to store # open parenthesis - # close parenthesis
                p_diff = 0;

                % go through each character of the remaining items to find when the parentheses resolve
                
                % go through each item
                for jj= ii+1:length(condition_text)
                    % get the text
                    item_text = condition_text{jj};
                    % crawl through, counting up the number of parenthesis
                    for kk=1:length(item_text)
                        if item_text(kk) == '('; p_diff = p_diff + 1;
                        elseif item_text(kk) == ')'; p_diff = p_diff - 1; end
                        % if the parenthesis resolve, stop searching
                        if p_diff == 0; break; end
                    end
                    if p_diff == 0; break; end
                    % otherwise, keep going
                end
                
                % if the loop ended without finding all parenthesis, give an error
                if p_diff ~=0
                    error('could not resolve parenthesis in index file.')
                end
                
                % return the value
                field_value = {condition_text{ii+1:jj}};
                
                
            else
                % if there are no special characters, just get the value of the field
                field_value = condition_text{ii+1};
            end
            
            

            % eliminate ", (, )
            field_value = strrep(field_value,'"','');
            field_value = strrep(field_value,'(','');
            field_value = strrep(field_value,')','');
            
            % save
            results.(desired_fields{dd}) = field_value;
            break
        end
    end
    % if it wasn't found, set the field to an empty string
    if ii == length(condition_text)
        results.(desired_fields{dd}) = '';
    end
end



function matlab_path = convert_path(obvius_path, fieldname)
% convert to matlab path format

% replace : with /
matlab_path = strrep(obvius_path,':','/');

% remove the ')'
matlab_path = strrep(matlab_path, ')', '');

% if it doesn't start with a slash, prepend server path
if ~strcmp(matlab_path(1),'/')
    if strcmp(fieldname, 'movie_xml_path')
        matlab_path = [moviexml_path matlab_path];
    end
end



function matlab_value = convert_value(obvius_value)
% convert to matlab format

% remove the ':'
matlab_value = strrep(obvius_value,':','');

% remove the ')'
matlab_value = strrep(matlab_value, ')', '');

% translate to number, if relevant
if ~isnan(str2double(matlab_value)); matlab_value = str2double(matlab_value);end




function cell_types = get_cell_types(condition_text)
% get cell types defined in condition text
% based on previous function load_cell_types.
%
% comments, authorship from previous code:
%
% Parses through Obvius index file to grab identified cell
% IDs. Assumes ':' identified new cell types. Does not recognize
% comments beginning with ';' well. Use at your own peril.
%
% shlens 2005-03-10 edits
% greschner 2006-02-10
% greschner 2008-09-23 search for condition check for equal:if isequal(textfile{i},condition)
%


i=1;
cell_types=struct;

while i <= length(condition_text)

    % found 'cell-types' keyword
    if ~isempty(strfind(condition_text{i},'cell-types'));
        i=i+1;

        % scan until find 'CONDITION' or 'condition'
        while (i<length(condition_text)) && ...
                (isempty(strfind(condition_text{i},'CONDITION'))) && ...
                (isempty(strfind(condition_text{i},'condition')))

            % scan until found keyword beginning with ":" (e.g. :on-y)
            if ~isempty(strfind(condition_text{i},':'));

                % format identified field name
                field_name = strrep(condition_text{i},':','');
                field_name = strrep(field_name,'(','');
                field_name = strrep(field_name,'-','_');

                % create blank structure field name
                cell_types.(field_name)=[];

                % scan until ')' or end of condition_text
                while i<=length(condition_text) && isempty(strfind(condition_text{i},')'))

                    % format identified field values (cell ID's)
                    field_value = strrep(condition_text{i},'(','');
                    field_value = strrep(field_value,')','');

                    % set identified field value
                    cell_types.(field_name)=[cell_types.(field_name) str2num(field_value)];
                    i=i+1;
                end

                % save last value
                % -- format identified field value (cell ID's)
                field_value = strrep(condition_text{i},'(','');
                field_value = strrep(field_value,')','');

                % -- set identified field value
                cell_types.(field_name)=[cell_types.(field_name) str2num(field_value)];

            end
            i=i+1;
        end

        i=length(condition_text); % stop loop
    end

    i=i+1;
end




% OLD CODE
%
%
% % set list of rrs fields (except rrs_prefix)
% rrs_fields = {'rrs_neurons_path','rrs_params_path','rrs_ei_path','rrs_sta_path'};
%
% % get rrs fields from index file
% index_paths = get_rrs_fields(condition_text, {'rrs_prefix',rrs_fields{:}});
% % index_paths.X = [] if not specified in index file
%
% % don't overwrite existing rrs_prefix, unless it's empty
% if ~isfield(datarun.names,'rrs_prefix') || isempty(datarun.names.rrs_prefix)
%     % load it from the index file (if none was specified, this will be empty)
%     datarun.names.rrs_prefix = index_paths.rrs_prefix;
% end
%
% % enter file names for rrs files, either as specified in index file, or generated from prefix
%
% % go through each field name
% for rr = 1:length(rrs_fields)
%     % don't overwrite existing field, unless it's empty
%     if ~isfield(datarun.names,rrs_fields{rr}) || isempty(datarun.names.(rrs_fields{rr}))
%         % load path from index file if it was specified
%         if ~isempty(index_paths.(rrs_fields{rr}))
%              datarun.names.(rrs_fields{rr}) = index_paths.(rrs_fields{rr});
%         else % otherwise generate name from rrs_prefix
%             % only if the prefix is not empty
%             if ~isempty(datarun.names.rrs_prefix)
%                 datarun.names.(rrs_fields{rr}) = [datarun.names.rrs_prefix '.' strrep(strrep(rrs_fields{rr},'_path',''),'rrs_','')];
%                     %strrep('rrs_','',strrep('_path','',rrs_fields{rr}))];
%             end
%         end
%     end
% end

