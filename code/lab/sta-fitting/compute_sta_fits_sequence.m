function datarun = compute_sta_fits_sequence(datarun, cell_specification, varargin)
% COMPUTE_STA_FITS_SEQUENCE     compute the fits for stas by sequentially 
%                               fitting temporal, center and surround.
%
% usage:  datarun = compute_sta_fits_sequence(datarun, cell_specification, varargin)
%
% INPUTS:
%   datarun                     A datarun structure that contains STAs
%   cell_specification          standard cell specification
%
% OPTIONAL INPUTS:
%     fit_instructions    []      A structure containing fit instructions,
%                                 see FIT_STA_SEQUENCE.M
%     verbose          false         
%   
%
%
% OUPUTS
%   datarun                     a datarun structure with fits results
%                               stored in datarun.matlab.sta_fits
%
%
%
% 2013-04 xyao
%



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
    datarun.matlab = [];
end

if ~isfield(datarun.matlab, 'sta_fits')
    temp_cell = cell(length(datarun.cell_ids), 1);
    datarun.matlab.sta_fits = temp_cell;
end


%% BODY OF THE FUNCTION

cell_indices = get_cell_indices(datarun, cell_specification);

num_rgcs = length(cell_indices);


% loop over cells and fit
for rgc = 1:num_rgcs
    
    if verbose
        fprintf('fitting the STA for cell %d... \n', datarun.cell_ids(cell_indices(rgc)))
    end
    temp_sta = datarun.stas.stas{cell_indices(rgc)};
    if isempty(fit_instructions)
        temp_fit = fit_sta_sequence(temp_sta);
    else
        temp_fit = fit_sta_sequence(temp_sta, fit_instructions);
    end
    
    datarun.matlab.sta_fits{cell_indices(rgc)} = temp_fit;

end
