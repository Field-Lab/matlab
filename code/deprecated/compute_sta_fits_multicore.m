function datarun = compute_sta_fits_multicore(datarun, cell_specification, varargin)
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
%   fit_params    []      A structure containing fit instructions,
%                               see FIT_STA.M
%
% OUPUTS
%   datarun                     a datarun structure with fits results
%                               stored in datarun.stas.stas.fits
%
%   COMPUTE_STA_FITS is a wrapper function for FIT_STA.M that conveniently
%   extracts STAs identified by cell_specification from a datarun 
%   structure, and fits these STAs according the fit_instructions
%   structure. See FIT_STA.M for details on how the STAs are fit.
%
% Author: GDF
% Date: 2011-06-11
%

% --- PARSE INPUTS ---
p = inputParser;

p.addRequired('datarun', @isstruct);
p.addRequired('cell_specification');

p.addParamValue('verbose', false, @islogical);
p.addParamValue('multicore_path', '~/matlab/multicore/', @ischar);
p.addParamValue('useWaitbar', true, @islogical);
p.addParamValue('fit_params', [], @isstruct);

p.parse(datarun, cell_specification, varargin{:});

verbose = p.Results.verbose;
settings.multicoreDir = p.Results.multicore_path;
settings.useWaitbar = p.Results.useWaitbar;
fit_params = p.Results.fit_params;

% check that stas have been loaded into datarun
if ~isfield(datarun.stas, 'stas')
    error('STAs are not contained in datarun, see LOAD_STA.M')
end


% --- INIITALIZE OUTPUTS ---

% check whether the fields 'matlab' and sta_fits' have been generated
if ~isfield(datarun, 'matlab')
    datarun = setfield(datarun, 'matlab', []);
elseif ~isfield(datarun.matlab, 'sta_fits')
    temp_cell = cell(length(datarun.cell_ids), 1);
    datarun.matlab.sta_fits = temp_cell;
end

if isempty(fit_params)
    fit_params.fit_surround = false;
end

%% --- BODY OF FUNCTION ----

% identify cells to use
cell_indices = get_cell_indices(datarun, cell_specification);
num_rgcs = length(cell_indices);

% put all stas to fit into a cell array
fit_info = cell(num_rgcs,1);
for rgc = 1:num_rgcs
    all_stas_to_fit{1} = datarun.stas.stas{cell_indices(rgc)};
    all_stas_to_fit{2} = fit_params;
    fit_info{rgc} = all_stas_to_fit;
end


% use multicore package to execute fits
fit_results_cell = startmulticoremaster(@fit_sta, fit_info, settings);



% put information into datarun.matlab.sta_fits    
for rgc = 1:num_rgcs 
    datarun.matlab.sta_fits{cell_indices(rgc)} = fit_results_cell{rgc};
end

        
