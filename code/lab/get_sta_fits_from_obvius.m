function datarun = get_sta_fits_from_obvius(datarun, cell_spec)
% get_sta_fits_from_obvius     Convert sta fits from obvius, and load in datarun.stas.fits
%
% usage:  datarun = get_sta_fits_from_obvius(datarun, cell_spec)
%
% arguments:  datarun - datarun struct with field datarun.obvius.sta_fits
%           cell_spec - which cells (see get_cell_indices for options)
%                       not specified - all
%
% outputs:    datarun - datarun struct with results stored in datarun.stas.fits
%
%
% 2009-04  gauthier
%



% ensure field exists
if ~isfield(datarun,'obvius')
    error('obvius info not loaded.  try load_obvius_sta_fits')
end
if ~isfield(datarun.obvius,'sta_fits')
    error('obvius STA fits not loaded.  use load_obvius_sta_fits')
end


if nargin==2
    % get cell indices
    cell_indices = get_cell_indices(datarun,cell_spec);

    % loop through cells
    for cc = 1:length(cell_indices)

        % get cell index and id
        cell_index = cell_indices(cc);

        % load converted fit
        datarun.stas.fits{cell_index} = sta_fit_from_obvius_fit(datarun.obvius.sta_fits{cell_index});
    end
else
    % loop through cells
    for cc = 1:length(datarun.obvius.sta_fits)
        % load converted fit
        datarun.stas.fits{cc} = sta_fit_from_obvius_fit(datarun.obvius.sta_fits(cc));
    end
end

