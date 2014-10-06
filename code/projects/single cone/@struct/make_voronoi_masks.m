function datarun = make_voronoi_masks(datarun, varargin)
% @STRUCT/MAKE_VORONOI_MASKS    Convert Voronoi cone tesselation to logical pixel masks
% usage: datarun = make_voronoi_masks(datarun, [opts])
%
% inputs: DATARUN
%
% opts: V,C             From running VORONOIN, if not given, VORONOIN will be run on
%                       DATARUN.CONES.CENTERS
%       WIDTH, HEIGHT   Dimensions of mask field.  If not given, taken from
%                       DATARUN.STIMULUS.FIELD_WIDTH and FIELD_HEIGHT
%       SCALE           Factor to scale the masks, for example if the cone
%                       finding was done on a 2x2 run but the voronoi
%                       should be a 1x1, SCALE would be 2
%
% outputs: Saves masks stack (HEIGHT x WIDTH x #CONES) to
%          DATARUN.CONES.MOSAIC.VORONOI_MASKS
%
% Takes polygon data from Voronoi tesselation of cone centers and converts 
% to a stack of logical pixel masks saved to DATARUN.
%
% See also: VORONOIN
%
% phli 2011-01
%
% 

opts = inputParser;
opts.addParamValue('V', []);
opts.addParamValue('C', []);
opts.addParamValue('scale', []);
opts.addParamValue('width', []);
opts.addParamValue('height', []);
opts.parse(varargin{:});
opts = opts.Results;

if isempty(opts.V) || isempty(opts.C)
    [opts.V,opts.C] = voronoin(datarun.cones.centers);
end

if isempty(opts.scale)
    if isfield(datarun.cones, 'mosaic') && isstruct(datarun.cones.mosaic) && ...
            isfield(datarun.cones.mosaic, 'voronoi_masks_scale') && ~isempty(datarun.cones.mosaic.voronoi_masks_scale)
        opts.scale = datarun.cones.mosaic.voronoi_masks_scale;
    elseif isfield(datarun, 'stimulus') && isfield(datarun.stimulus, 'stixel_width') && ...
            isfield(datarun.stimulus, 'stixel_height') && datarun.stimulus.stixel_width == datarun.stimulus.stixel_height
        opts.scale = datarun.stimulus.stixel_width;
    else
        opts.scale = 1;
    end
end

if isempty(opts.width)
    opts.width = datarun.stimulus.field_width .* opts.scale;
end

if isempty(opts.height)
    opts.height = datarun.stimulus.field_height .* opts.scale;
end

datarun.cones.mosaic.voronoi_masks = make_voronoi_masks(opts.V, opts.C, opts.width, opts.height, opts.scale);
datarun.cones.mosaic.voronoi_masks_scale = opts.scale;