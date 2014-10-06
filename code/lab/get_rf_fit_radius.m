function radius = get_rf_fit_radius(datarun, cell_spec, varargin)
% get_rf_fit_radius     This function calculated the rf fit radius of cells in 
%                       in datarun specified by cell_spec
%
% usage:  radius = get_rf_fit_radius(datarun, cell_spec, varargin)
%
% arguments:     datarun - datarun structure
%            cell_spec - specifies cells: see get_cell_indices
%
% outputs:     radius - radii for the cells specified by cell_spec
%
%
% optional params, their default values, and what they specify:
%
% fits_to_use           vision          'vision' or 'obivus' fits can be specified
% units                 pixels          'pixels of 'microns' can be specified
% microns_per_pixel       5.5            For results to be trusted absolutely, this value needs to be
%                                        set according to the datarun
%
% 2010-04 GDF
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('fits_to_use', 'vision');
p.addParamValue('units', 'pixels');
p.addParamValue('microns_per_pixel', 5.5);
p.parse(varargin{:});

fits_to_use = p.Results.fits_to_use;

% check that the datarun.stimulus field is filed
if ~isfield(datarun, 'stimulus')
    error('Stimulus information is not contained in datarun: load stimulus information');
end

% check that stixel height and width are the same
if datarun.stimulus.stixel_height ~= datarun.stimulus.stixel_width
    error('Stixel height and width are not equal: function assumes square pixels and uses stixel height')
end

% BODY OF FUNCTION

% indices of cells in cell_spec
cell_indices = get_cell_indices(datarun, cell_spec);

% number of cells to calculate rf radius
num_cells = length(cell_indices);

% intialize output
radius = zeros(num_cells,1);


for cc = 1:num_cells

    switch fits_to_use
    
        case 'vision'
            sta_fit = datarun.vision.sta_fits{cell_indices(cc)};
        case 'obvius'
            sta_fit = datarun.obvius.sta_fits{cell_indices(cc)};
        otherwise
            error('Unknown fits_to_use: must be vision or obvius')

    end

    sd = sta_fit.sd;
    radius(cc) = geomean(sd) * datarun.stimulus.stixel_height;
    
    % if specified, put radius in units of microns
    if strcmp(p.Results.units, 'microns')
        radius(cc) = radius(cc) * p.Results.microns_per_pixel;
    end

end

    
