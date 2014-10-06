function datarun = load_ei(datarun, cell_specification, varargin)
% LOAD_EI     Load information from an ei file
%
% usage:  datarun = load_ei(datarun, cell_specification, params)
%
% arguments:  datarun - datarun struct with field specifying the neurons file
%                         (datarun.names.rrs_ei_path - absolute path to ei file,
%                          e.g. '/Analysis/Greschner/2005-04-26-1/data000/data000.ei')
%
%              cell_specification: Vision ID of cell, not matlab ID
%
%              params - array_type = 61 | 512 | 519  default 512, uses fake id
%                       array_id  default [], overrides array_type
%                       zero_disconnected, default false
%                       keep_java_ei, default true for backward
%                           compatibility but better to use false
%                       require  if true, ei is NOT loaded if it's already stored in datarun
%                                if false, ei is loaded no matter what
%
%
% outputs:    datarun - datarun struct with the following fields added, as possible
%
%     	datarun.ei.eis        - N-length cell array, each entry an ExF matrix (E = # electrodes, F = # frames)
%       datarun.ei.positions  - Ex2 matrix of x-y coordinates of electrode positions
%       datarun.ei.java_ei    - java ei file object (if keep_java_ei param is true)
%
% greschner changed from shlens 2006-07-23
% tamachado 2009-06-15
% gauthier  2009-10-23, cleaned up, added storage of java_ei
% gauthier  2010-02, array info now obtained from load_array_info
% phli      2011-08, now loads number of spikes used for each EI
% calculation
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addParamValue('array_type', 512, @(x) any(x == [  61 512 519  ]));%
p.addParamValue('array_id', [], @(x)x>0 && mod(x,1)==0)
p.addParamValue('zero_disconnected', 0);%
p.addParamValue('keep_java_ei', true);
p.addParamValue('require', true);%
% parse inputs
p.parse(varargin{:});
params = p.Results;


% if no array_id is specified...
if isempty(params.array_id)
    
    % look in datarun for an array_id
    if isfield(datarun,'piece') && isfield(datarun.piece,'array_id') && ~isempty(datarun.piece.array_id)
        params.array_id = datarun.piece.array_id;
    else
        % Try globals file
        datarun = load_globals(datarun);
        if isfield(datarun, 'globals') && datarun.globals.getImageCalibrationParams.arrayID
            params.array_id = datarun.globals.getImageCalibrationParams.arrayID;
        else
            % if datarun.piece.array_id was empty, choose a valid array_id using array_type
            params.array_id = fake_array_id(params.array_type);
            fprintf('!!!!!! ARRAY ID WAS NOT FOUND IN DATARUN.  IT WAS AUTOMATICALLY SET TO %d !!!!!!',params.array_id)            
        end
        
        % and save this value in datarun
        datarun.piece.array_id = params.array_id;
        
    end
else
    % if an array id is specified, ensure it agrees with datarun
    if isfield(datarun,'piece') && isfield(datarun.piece,'array_id')
        if datarun.piece.array_id ~= params.array_id
            error('the specified array id (%d) does not match the array id from datarun (%d)',params.array_id,datarun.piece.array_id)
        end
    else
        % if datarun.piece.array_id is empty, use the specified array id
        datarun.piece.array_id = params.array_id;
        fprintf('!!!!!! ARRAY ID WAS NOT FOUND IN DATARUN.  IT WAS AUTOMATICALLY SET TO %d !!!!!!',params.array_id)
    end
end



% but if array_id IS specified, identify the corresponding array_type
% if params.array_id < 500
%     params.array_type = 61;
% elseif params.array_id >= 500 && params.array_id < 1500
%     params.array_type = 512;
% elseif params.array_id >= 1500 && params.array_id < 2500
%     params.array_type = 519;
% else
%     error('array_id %d not recognized',params.array_id)
% end





% ensure correct field exists

% if the field is not there...
if ~isfield(datarun,'ei') || ~isfield(datarun.ei,'eis')
    % initialize it as a cell array
    datarun.ei.eis = cell(length(datarun.cell_ids),1);
else
    % or if it is too small...
    if length(datarun.ei.eis) < length(datarun.cell_ids)
        % ... make it long enough
        datarun.ei.eis{length(datarun.cell_ids)} = [];
    end
end






% IDENTIFY EIs TO LOAD


% get cell indices
cell_indices = get_cell_indices(datarun,cell_specification);

% determine whether any EIs need to be loaded
if params.require && isfield(datarun,'ei') && isfield(datarun.ei,'eis')
    % if require is true, only load up EIs which are not already loaded
    to_be_loaded = [];
    for cc=1:length(cell_indices)
        if isempty(datarun.ei.eis{cell_indices(cc)})
            to_be_loaded = [to_be_loaded cc];
        end
    end
else
    to_be_loaded = 1:length(cell_indices);
end

% if there are no EIs to be loaded, abandon ship
%if isempty(to_be_loaded)
%    return; end






% LOAD EI FILE AND ARRAY INFO

% check directory
if ~exist(datarun.names.rrs_ei_path,'file')
    error('load_ei: bad directory path: ''%s''',datarun.names.rrs_ei_path);
end

% load ei object into datarun if not already present
if ~isfield(datarun,'ei') || ~isfield(datarun.ei,'java_ei') || isempty(datarun.ei.java_ei)
    datarun.ei.java_ei = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(datarun.names.rrs_ei_path);
end

% load array info
array_info = load_array_info(datarun,2);

% create java electrode map
%electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(params.array_id);






% LOAD OTHER INFORMATION

% note number of right and left points
datarun.ei.nrPoints = datarun.ei.java_ei.nrPoints;
datarun.ei.nlPoints = datarun.ei.java_ei.nlPoints;

% note electrode positions
datarun.ei.position = array_info.positions;
%datarun.ei.position = electrode_positions(params.array_type);

% note array boundary, useful for setting axis limits when plotting
datarun.ei.array_bounds_x = array_info.x_bounds;
datarun.ei.array_bounds_y = array_info.y_bounds;

datarun.ei.disconnected = array_info.disconnected;




% Note number of spikes used for each cell and the average number of spikes
datarun.ei.num_spikes = zeros(size(datarun.cell_ids));
for i = 1:numel(datarun.cell_ids)
    datarun.ei.num_spikes(i) = datarun.ei.java_ei.getNSpikes(datarun.cell_ids(i));
end



% LOAD ONE EI AT A TIME

% loop through each cell
for cc = to_be_loaded

    % get cell index
    cell_index = cell_indices(cc);

    % load EI data
    ei = datarun.ei.java_ei.getImage(datarun.cell_ids(cell_index));

    % throw out the error
    ei = squeeze(ei(1,:,:));

    % throw out the triggers
    ei = ei(2:end,:);

    % grab paramters
    electrodes = size(ei,1);
    frames = size(ei,2);

    % zero disconnected electrodes
    if params.zero_disconnected
        %         for i=1:electrodes
        %             if electrodeMap.isDisconnected(i)==1
        %                 ei(i,:) = 0;
        %             end
        %         end
        ei(datarun.ei.disconnected,:) = 0;
    end

    % restore double
    datarun.ei.eis{cell_index} = double(ei);

end



if ~params.keep_java_ei
    datarun.ei.java_ei.close();
    datarun.ei = rmfield(datarun.ei, 'java_ei');
end



% switch params.array_type
% 
%     case 61
%         datarun.ei.array_bounds_x = 350*[-1 1];
%         datarun.ei.array_bounds_y = 300*[-1 1];
% 
%     case 512
%         datarun.ei.array_bounds_x = 1100*[-1 1];
%         datarun.ei.array_bounds_y = 500*[-1 1];
% 
%     case 519
%         datarun.ei.array_bounds_x = 900*[-1 1];
%         datarun.ei.array_bounds_y = 800*[-1 1];
% 
%     otherwise
%         error('bug in code? params.array_type should have been error checked above')
% end





% old code to get electrode positions

% % import java classes
% import('edu.ucsc.neurobiology.vision.io.*');
% import('edu.ucsc.neurobiology.vision.electrodemap.*');
%
% % create java electrode map
% electrodeMap = ElectrodeMapFactory.getElectrodeMap(params.array_id);
% 
% % grab positions
% for i=1:(datarun.ei.java_ei.nElectrodes-1)
%     position(i,1) = electrodeMap.getXPosition(i);
%     position(i,2) = electrodeMap.getYPosition(i);
% end
% datarun.ei.position = double(position);

