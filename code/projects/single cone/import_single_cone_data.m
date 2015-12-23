function varargout = import_single_cone_data(datarun,path,varargin)
% import_single_cone_data     read text files with information about single cone sampling
%
% usage:  varargout = import_single_cone_data(datarun,path)
%
% arguments:    datarun - datarun struct
%                  path - path to text files.  if starts with '/', will be taken as absolute.
%                           otherwise, will be assumed relative to this path:
%                           '/snle/lab/Experiments/Array/Shared/one/'
%
% optional arguments
%               overwrite_cell_types  false     overwrite cell type 
%               load_raw_results      false     loads raw results struct (mostly for getting likelihood surface 
%                                               from Bayesian analyses)
%
% outputs:
%
%
%  if datarun is empty, output is a struct with the following fields
%
%       cell_ids            list of RGC ids
%       cell_types          cell type of each RGC
%           (0 = none, 1 = on-parasol, 2 = off-parasol, 3 = on-midget, 4 = off-midget, 5 = sbc, 6+ = other)
%       rgc_COMs            center of mass of each RGC RF
%       cells_in_roi        which RGCs are in the region of interest
%       rgc_fit_info        Nx8 matrix, parameters of DOG fit to each RGC RF
%           1,2 = center location, 3=center radius, 4=center amplitude, 5=surround radius, 6=surround amplitude,
%           7=radius of cones used in the fit (units=stixels), 8=error of the fit
%
%       cone_centers        Nx2 matrix of x-y cone center locations
%       cone_ids            list of cone IDs
%       cone_rgb_values     Nx3 matrix of the r-g-b sensitivity of each cone
%       cone_types          string vector giving type of each cone (L,M,S,U for unknown)
%       cone_roi            which cones are in the region of interest
%       cone_likelihoods    likelihood ratio (from EM cone type classification)
%
%       cone_weights        weight matrix of connection strengths (rows = cones, columns = cells)
%
%
%
%
%  if datarun is a struct, output is the same struct with data entered in these fields
%
%   datarun.cones.centers   - Nx2 matrix with x,y coordinates
%                .types     - N-length char vector, 'L','M','S', or 'U' for unknown
%                .rgb       - Nx3 matrix withrgb triplet for each cone
%                .weights   - NxM matrix, each column gives the cone weights for one RGC
%                .rf_fits   - Mx1 cell array of structs, each with fields
%                               center
%                               center_radius
%                               center_scale
%                               surround_radius
%                               surround_scale
%                               fit_radius
%                               error
%
%
%
% examples:
%
%   [datarun,extras] = import_single_cone_data(datarun);
%
%   datarun = import_single_cone_data(datarun,'kiwi');
%
%   cone_data = import_single_cone_data([],'kiwi');
%
%   cone_data = import_single_cone_data([],'/snle/home/gauthier2/Desktop/cone_data/');
%
%
%
% 2008-10 gauthier
% 2009-01 gauthier, added support for cone classification fields
% 20??-?? phli, added optional loading of raw results (mostly for getting likelihood surface from Bayesian analyses)
%

p = inputParser;
p.addParamValue('overwrite_cell_types', false, @islogical);
p.addParamValue('load_raw_results',     false, @islogical);
p.parse(varargin{:});

% get file prefix

% if no path specified, use the nickname or short name as the specification
if nargin == 1
    if isfield(datarun,'names')
        if isfield(datarun.names,'nickname')
            path = datarun.names.nickname;
        elseif isfield(datarun.names,'short_name')
            path = datarun.names.short_name;
        end
    end
end

if ~exist('path','var')
    error('Path to text files was not specified, and datarun.names did not have fields ''nickname'' or ''short_name''.')
end
 

% if path is relative
if ~strcmp(path(1),'/') && ~strcmp(path(1), '~');
    
    path2data=datarun.names.rrs_prefix;
    tmp=regexp(path2data,'data');
    path2data=path2data(1:tmp(1)-1);
    
    % assume it is in the expected folder
    file_prefix = [path2data path '/'];
%     file_prefix = [single_cone_path path '/'];
else
    % otherwise, use the full path
    file_prefix = [path '/']; % add the trailing slash, just to be kind
end



% load files

rgcs_file = dlmread([file_prefix 'rgcs.txt'],'\t');
cones_file = dlmread([file_prefix 'cones.txt'],'\t');
connections_file = dlmread([file_prefix 'connections.txt'],'\t');

fprintf('Reading text files from %s:\n',file_prefix)
fprintf('rgcs:        [%d x %d] matrix.\n',size(rgcs_file,1),size(rgcs_file,2))
fprintf('cones:       [%d x %d] matrix.\n',size(cones_file,1),size(cones_file,2))
fprintf('connections: [%d x %d] matrix.\n',size(connections_file,1),size(connections_file,2))

if p.Results.load_raw_results
    results = load([file_prefix 'results.mat']);
end


% put file data into variables


% rgc info

% id number for each RGC
cell_ids = rgcs_file(:,1);

% which cell type (0 = none, 1 = on-parasol, 2 = off-parasol, 3 = on-midget, 4 = off-midget, 5 = sbc, 6+ = other)
cell_types = rgcs_file(:,2);

% whether a cell is in the region of interest (1 = in ROI, 0 = not in ROI)
cells_in_roi = rgcs_file(:,3);

% RGC center, based on center of mass of the cones
% first column is x coordinate, second column is y coordinate
rgc_COMs = rgcs_file(:,[4 5]);

% DOG fit
rgc_fit_info = rgcs_file(:,[6 7 8 9 10 11 12 13]);



% cone info

% id number for each cone
cone_ids = cones_file(:,1);

% center location of each cone
cone_centers = cones_file(:,[2 3]);

% color of each cone (L,M,S, or U for unknown)
cone_types = cone_type_nums_to_char( cones_file(:,4) );

% measured red-green-blue color sensitivity for each cone
cone_rgb_values = cones_file(:,[5 6 7]);

% cone ROI
cone_roi = cones_file(:,8);

% identity derived from EM cone classification
cone_types_em = cone_type_nums_to_char( cones_file(:,9) );

% likelihood ratio for each cone
cone_likelihoods = cones_file(:,10);

% identify derived from kmeans cone classification
cone_types_kmeans = cone_type_nums_to_char( cones_file(:,11) );


% weights info

% C by R matrix (C = # of cones, R = # of ganglion cells), giving cone weights for each RGC
cone_weights = connections_file;




% return output
if isstruct(datarun)
    % load data into datarun

    datarun.cones.centers = cone_centers;
    datarun.cones.types = cone_types;
    datarun.cones.rgb = cone_rgb_values;
    datarun.cones.likelihoods = cone_likelihoods;
    datarun.cones.types_em = cone_types_em;
    datarun.cones.types_kmeans = cone_types_kmeans;
    datarun.cones.roi = cone_roi;

    % get cell indices of the RGCs in these data
    cell_indices = get_cell_indices(datarun,cell_ids);
    
    % initialize weight matrix with zeros
    datarun.cones.weights = zeros(length(cone_ids),length(datarun.cell_ids));
    % then fill in entries
    datarun.cones.weights(:,cell_indices) = cone_weights;

    % enter fits and center points
    for cc = 1:length(cell_indices)
        clear the_fit
        % if it's all zeros, then no fit was made, so set it to empty
        if all(rgc_fit_info(cc,:)==0)
            the_fit = [];
        else
            % otherwise, enter all the numbers
            the_fit.center = rgc_fit_info(cc,1:2);
            the_fit.center_radius = rgc_fit_info(cc,3);
            the_fit.center_scale = rgc_fit_info(cc,4);
            the_fit.surround_radius = rgc_fit_info(cc,5);
            the_fit.surround_scale = rgc_fit_info(cc,6);
            the_fit.fit_radius = rgc_fit_info(cc,7);
            the_fit.error = rgc_fit_info(cc,8);
        end
        
        % get the cell_index of this cell
        cell_index = cell_indices(cc);
        cell_index = get_cell_indices(datarun,datarun.cell_ids(cell_indices(cc)));
        
        % save this fit in datarun
        datarun.cones.rf_fits{cell_index} = the_fit;
     
        
        % if the center point is empty
        if isfield(datarun,'stas') && isfield(datarun.stas, 'rf_coms')
            if isempty(datarun.stas.rf_coms{cell_index})
                % load the center point from the text file
                datarun.stas.rf_coms{cell_index} = rgc_COMs(cc,:);
            end
        end
    end

    if p.Results.overwrite_cell_types
        for c_type = 1:6
            temp_inds = find(cell_types == c_type);
            datarun.cell_types{c_type}.cell_ids = cell_ids(temp_inds);
        end
        warning('overwriting cell type information')
    end

    
    
    % return datarun
    varargout(1) = {datarun};
    
    % if desired, also return extras
    if nargout > 1
        extras.cell_ids = cell_ids;
        extras.rgc_roi = cells_in_roi;
        
        if p.Results.load_raw_results
            extras.results = results;
        end
        
        varargout(2) = {extras};
    end
     

elseif isempty(datarun)

    % return variables in a struct

    % clean up workspace
    clear rgcs_file cones_file connections_file p datarun path file_prefix varargin

    % get list of variables
    var_list=whos;

    % put them in a struct
    output = struct;
    for vv = 1:length(var_list)
        eval(['output.' var_list(vv).name ' = ' var_list(vv).name ';']);
    end
    
    % return the struct
    varargout(1) = {output};

else
    error('second argument must either be empty, or a datarun struct.')
end


function chars = cone_type_nums_to_char(types)

chars = char;
chars(types == 1,1) = 'L';
chars(types == 2,1) = 'M';
chars(types == 3,1) = 'S';
chars(types == 0,1) = 'U';
