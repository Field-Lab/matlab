function weights = dog_fit(locations, varargin)

% provide a list of locations (for cones) and a set of fit parameters in
% varargin, and the function returns a set of weights from the DOG fit to
% the cell.


p = inputParser;

p.addRequired('locations');
p.addParamValue('center', [0 0]);
p.addParamValue('center_radius', 1);
p.addParamValue('center_scale', 1);
p.addParamValue('surround_radius', 1);
p.addParamValue('surround_scale', 0);
p.addParamValue('fit_radius', 100);
p.addParamValue('error', 0);

p.parse(locations, varargin{:});

center = p.Results.center;
center_radius = p.Results.center_radius;
center_scale = p.Results.center_scale;
surround_radius = p.Results.surround_radius;
surround_scale = p.Results.surround_scale;

% this function assumes the fit is circular
distances = ipdm(locations, center);

if center_radius ~= 0
    center_Gauss_weights = center_scale .* exp(((distances/center_radius).^2) .* -0.5);
end
if surround_radius ~= 0
    surround_Gauss_weights = surround_scale .* exp(((distances/surround_radius).^2) .* -0.5);
else
    surround_Gauss_weights = zeros(length(distances),1);
end
weights = center_Gauss_weights - surround_Gauss_weights;







