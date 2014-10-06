function datarun = get_sta_fits_from_vision(datarun, cell_spec)
% get_sta_fits_from_vision     Convert sta fits from vision, and load in datarun.stas.fits
%
% usage:  datarun = get_sta_fits_from_vision(datarun, cell_spec)
%
% arguments:  datarun - datarun struct with field datarun.vision.sta_fits
%           cell_spec - which cells (see get_cell_indices for options)
%                           if not specified, set to 'all'
%
% outputs:    datarun - datarun struct with results stored in datarun.stas.fits
%
%
% 2009-04  gauthier
%



% ensure field exists
if ~isfield(datarun,'vision')
    error('Vision info not loaded yet')
end
if ~isfield(datarun.vision,'sta_fits')
    error('Vision STA fits not loaded yet')
end

% by default, load all
if ~exist('cell_spec','var')
    cell_spec = 'all';
end


% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);

% loop through cells
for cc = 1:length(cell_indices)
    
    % get cell index and id
    cell_index = cell_indices(cc);
    
    % load converted fit
    datarun.stas.fits{cell_index} = sta_fit_from_vision_fit(datarun,datarun.vision.sta_fits{cell_index});
end

