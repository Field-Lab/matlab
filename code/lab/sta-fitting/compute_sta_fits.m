function datarun = compute_sta_fits(datarun, cell_specification, varargin)
%
% COMPUTE_STA_FITS  computes the fits for the stas of cells specified by
%                    cell_specification
%
%   datarun = compute_sta_fits(datarun, cell_specification, varargin)
%
% INPUTS:
%   datarun                     A datarun structure that contains STAs
%   cell_specification          standard cell specification
%   
% OPTIONAL INPUTS:
%   fit_instructions    []      A structure containing fit instructions,
%                               see FIT_STA.M
%
% OUPUTS
%   datarun                     a datarun structure with fits results
%                               stored in datarun.matlab.sta_fits
%
%   COMPUTE_STA_FITS is a wrapper function for FIT_STA.M that conveniently
%   extracts STAs identified by cell_specification from a datarun 
%   structure, and fits these STAs according the fit_instructions
%   structure. See FIT_STA.M for details on how the STAs are fit.
%
% Author: GDF
% Date: 2011-06-11

p = inputParser;

p.addRequired('datarun', @isstruct);
p.addRequired('cell_specification');

p.addParamValue('fit_instructions', [], @isstruct);
p.addParamValue('verbose', false, @islogical);

p.parse(datarun, cell_specification, varargin{:});

fit_instructions = p.Results.fit_instructions;
verbose = p.Results.verbose;

% check that stas have been loaded into datarun
if ~isfield(datarun.stas, 'stas')
    error('STAs are not contained in datarun, see LOAD_STA.M')
end

% initialize_output
if ~isfield(datarun, 'matlab')
    datarun = setfield(datarun, 'matlab', []);
elseif ~isfield(datarun.matlab, 'sta_fits')
    temp_cell = cell(length(datarun.cell_ids), 1);
    datarun.matlab.sta_fits = temp_cell;
end


%% BODY OF FUNCTION
cell_indices = get_cell_indices(datarun, cell_specification);

num_rgcs = length(cell_indices);

% loop over cells and fit
for rgc = 1:num_rgcs
    
    if verbose
        fprintf('fitting the STA for cell %d... \n', datarun.cell_ids(cell_indices(rgc)))
    end
    
    temp_sta = datarun.stas.stas{cell_indices(rgc)};
    if isempty(fit_instructions)
        temp_fit_params = fit_sta(temp_sta);
    else
        temp_fit_params = fit_sta(temp_sta, fit_instructions);
    end
    
    if isempty(temp_fit_params)
        temp_id = datarun.cell_ids(cell_indices(rgc));
        warn_message = ['cell ',num2str(temp_id), ' has no sig stixels and no fit'];
        warning(warn_message)
    end
    
    datarun.matlab.sta_fits{cell_indices(rgc)} = temp_fit_params;
    
end

        
