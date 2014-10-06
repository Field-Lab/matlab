function [crs crsx rasterhistx rasterhistsy] = calc_allcones_cr(datarun, cellids, triggers, varargin)
% CALC_ALLCONES_CR
% For allcones analysis based on contrast/response scale matching rather
% than simply comparing rates

opts = inputParser();
opts.addParamValue('template', []);
opts.addParamValue('baselines', []);
opts.addParamValue('boxstart', 0.05);
opts.addParamValue('boxend',   0.25);
opts.parse(varargin{:});
opts = opts.Results;

% Defaults to boxcar template
if isempty(opts.template)
    template.x = -0.1:0.02:0.5;
    template.y = zeros(size(template.x));
    template.y(template.x > opts.boxstart & template.x < opts.boxend) = 1;
    template.y = template.y ./ sum(template.y);
else
    template = opts.template;
end

stimulus = datarun.stimulus;
stimulus = parse_cr_rgbs(stimulus);
[crs rasterhistsy] = single_map_contrast_response(datarun, triggers, cellids, template, cell2mat(stimulus.cr.single_cones'));

% Process into nicer dimensions
crs = permute(reshape(crs, length(datarun.rgcs), [], stimulus.numcones), [3 2 1]);
rasterhistsy = permute(reshape(rasterhistsy, length(datarun.rgcs), [], stimulus.numcones), [3 2 1]);
crsx = cell2mat(stimulus.cr.single_cone_intensities);
rasterhistx = template.x;