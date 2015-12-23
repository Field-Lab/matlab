function blurred_cone_mosaic = compute_blurred_cone_mosaic(datarun, varargin)


% PARSE INPUTS
p = inputParser;

%p.addRequired('datarun', @isstruct)
p.addParamValue('blur_factor', 5, @isnumeric);
p.addParamValue('blur_amplitude', 1, @isnumeric);
p.addParamValue('saturation_point', 1.0, @isnumeric);
p.addParamValue('verbose', false', @islogical);

p.parse(varargin{:});


% -- BODY OF FUNCTION --

% extract paramters about field size and cone number
height = datarun.stimulus.field_height;
width = datarun.stimulus.field_width;
num_cones = length(datarun.cones.types);

% set the paramaters for generating the Gaussian blur function
blur_params.center_radius = p.Results.blur_factor;
blur_params.center_scale = p.Results.blur_amplitude;
blur_params.surround_radius = 0;
blur_params.surround_scale = 0;
blur_params.x_size = width;
blur_params.y_size = height;

% get indices to L and M cones (skip S)
L_indices = find(datarun.cones.types == 'L');
M_indices = find(datarun.cones.types == 'M');

% initialize some matrices
cone_mosaic = zeros(height,width,3);
colored_L = zeros(height,width,2);
colored_M = zeros(height,width,1);

if p.Results.verbose
    h = waitbar(0, 'looping over cones...');
end

% loop over cones, blur and images for each cone.
for cn = 1:num_cones
    
    blur_params.center = datarun.cones.centers(cn,:);
    blurred_cone = make_gaussian(blur_params);
    
    if ~isempty(find(L_indices == cn, 1))
        colored_blurred_cone = cat(3, blurred_cone, colored_L);
        cone_mosaic = cone_mosaic + colored_blurred_cone;
    elseif ~isempty(find(M_indices == cn, 1))
        colored_blurred_cone = cat(3, colored_M, blurred_cone, colored_M);
        cone_mosaic = cone_mosaic + colored_blurred_cone;
    end

    if p.Results.verbose
        waitbar(cn ./ num_cones)
    end
    
end    

close(h)

extremum_val = max(max(max(cone_mosaic)));
normalization_val = extremum_val ./ p.Results.saturation_point;
blurred_cone_mosaic = cone_mosaic ./ normalization_val;

blurred_cone_mosaic(blurred_cone_mosaic > 1.0) = 1.0;



