function datarun = load_params(datarun, varargin)
% LOAD_PARAMS     Load information from a vision params file
%
% usage:  datarun = load_params(datarun, <params>)
%
% arguments:  datarun - datarun struct with field specifying the params file
%                         (datarun.names.rrs_params_path - absolute path to params file,
%                          e.g. '/Analysis/Greschner/2005-04-26-1/data000/data000.params')
%              params - struct of optional parameters (see below)
%
% outputs:    datarun - datarun struct with the following fields added, as desired
%
% 	piece information:
%     	datarun.vision.cell_types
%     	datarun.vision.sta_fits
%
%
% optional fields in params, their default values, and what they specify:
%
% verbose            	false	print information about loading
% sync_cell_ids         true    if datarun.cell_ids doesn't match the params file, give an error
%                               if datarun.cell_ids doesn't exist, set it to match the params file
%
% load_cell_types       true    load cell types?
% cell_type_depth       2       how many levels deep to look for cell types in the classification ['all']
% order_cell_types      true    put cell types in the canonical order
%
% load_sta_fits         'all'  	which STA fits to load from the params file
%                                   loaded into datarun.vision.sta_fits
%                                   to load none, make []
%
% 2008     gauthier
% 2009?    greschner?  ?
% 2009-09  gauthier        added option to load STA fits
% 2010     greschner       increase cell_type_depth to 3 
% 2012-01  jepson, phli    close params file at the end, unless requested to keep
% 2012-05  phli            option to load Vision timecourses
% 2012-10  phli            loads fits automatically if possible
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('sync_cell_ids', true);
p.addParamValue('load_cell_types', true);
p.addParamValue('cell_type_depth', 3);
p.addParamValue('order_cell_types', true);
p.addParamValue('load_sta_fits', 'all');
p.addParamValue('load_timecourses', 'all');
p.addParamValue('keep_java_params', false);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% LOAD THE PARAMS FILE

% show what's going on
if params.verbose;
    fprintf('\nLoading Java paramsFile object...');
    start_loading = clock; % note when it started
end
% load rrs params file
java_params = edu.ucsc.neurobiology.vision.io.ParametersFile(datarun.names.rrs_params_path);
% display how long it took 
if params.verbose;
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_loading));
end



% GET CELL IDS FROM THE PARAMS FILE

% cell IDs
cell_ids_from_params = java_params.getIDList';

% if datarun.cell_ids doesn't match the STAs, give an error.
% if datarun.cell_ids doesn't exist, define it to match the STAs
if params.sync_cell_ids
    datarun = sync_cell_ids(datarun, cell_ids_from_params,...
        sprintf('params file %s',datarun.names.rrs_params_path));
end




% LOAD CELL TYPES

if params.load_cell_types

    % show what's going on
    if params.verbose;
        fprintf('\nLoading cell types...');
        start_loading = clock; % note when it started
    end

    % get cell types from java object
    cell_types = load_rrs_cell_types(java_params, struct('cell_type_depth',params.cell_type_depth));

    % order cell types, if desired
    if params.order_cell_types
        cell_types = order_cell_types(cell_types);
    end

    datarun = load_celltypes_into_datarun(datarun, cell_types, params.verbose);
    
    % display how long it took
    if params.verbose, fprintf(' done (%0.1f seconds)\n',etime(clock,start_loading)); end
end


% LOAD STA FITS (if desired)
if ~isempty(params.load_sta_fits)
    
    % note cell indices to load
    cell_indices = get_cell_indices(datarun,params.load_sta_fits);
    
    % ensure storage field exists...
    if ~isfield(datarun,'vision') || ~isfield(datarun.vision,'sta_fits')
        datarun.vision.sta_fits = cell(length(datarun.cell_ids),1);
    else
        % and has enough elements
        if length(datarun.vision.sta_fits) < length(datarun.cell_ids)
            datarun.vision.sta_fits{length(datarun.cell_ids)} = [];
        end
    end
    
    % note when it started
    if params.verbose;
        fprintf('Loading %d STA fits...',length(cell_indices))
        start_loading = clock;
    end
    
    % load fit for each cell
    for cc=1:length(cell_indices)
        
        % get cell index and cell id
        cell_index = cell_indices(cc);
        cell_id = datarun.cell_ids(cell_index);
        
        % initialize fit
        fit = struct;
        
        % load up each part of the fit
        try
            fit.mean(1)=getDoubleCell(java_params, cell_id, 'x0');
            fit.mean(2)=getDoubleCell(java_params, cell_id, 'y0');
            fit.sd(1)=getDoubleCell(java_params, cell_id, 'SigmaX');
            fit.sd(2)=getDoubleCell(java_params, cell_id, 'SigmaY');
            fit.angle=getDoubleCell(java_params, cell_id, 'Theta');

            % if fit did not converge, enter appropriate data
            if isnan(fit.angle)
                fit.mean = [-1 -1];
                fit.sd = [NaN NaN];
                fit.angle = NaN;
            end   
        catch
        end

        % store in datarun
        datarun.vision.sta_fits{cell_index}=fit;

    end
    
    % If no sta fits loaded into main location yet, go ahead and load these
    % No point if stas not loaded yet though, as load_sta will overwrite
    if isfield(datarun, 'stas') && ~isfield(datarun.stas, 'fits')
        datarun = get_sta_fits_from_vision(datarun);
    end
    
    % show how long it took
    if params.verbose, fprintf(' done (%0.1f seconds)\n',etime(clock,start_loading)); end
    
end


% Load timecourses if desired
if ~isempty(params.load_timecourses)
    cell_indices = get_cell_indices(datarun,params.load_sta_fits);

    % ensure storage field exists...
    if ~isfield(datarun,'vision') || ~isfield(datarun.vision,'timecourses')
        timecourses = struct([]);
    end

    for cc = 1:length(cell_indices)
        cell_index = cell_indices(cc);
        cid = datarun.cell_ids(cell_index);

        if java_params.hasParameter('RedTimeCourse')
            timecourses(cell_index).r = java_params.getArrayCell(cid, 'RedTimeCourse'); end
        if java_params.hasParameter('GreenTimeCourse')
            timecourses(cell_index).g = java_params.getArrayCell(cid, 'GreenTimeCourse'); end
        if java_params.hasParameter('BlueTimeCourse')
            timecourses(cell_index).b = java_params.getArrayCell(cid, 'BlueTimeCourse'); end

        if ~isempty(timecourses)
            datarun.vision.timecourses = timecourses; end
    end
end


if params.keep_java_params
    datarun.vision.java_params = java_params;
else
    % .close requires boolean argument specifying whether to save out the
    % params file before closing it
    java_params.close(false)
end