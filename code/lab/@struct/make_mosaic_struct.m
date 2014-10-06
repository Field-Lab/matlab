function datarun = make_mosaic_struct(datarun, varargin)
% @struct/make_mosaic_struct     datarun wrapper for make_mosaic_struct, respects existing fields
%
% usage:  datarun = make_mosaic_struct(datarun, varargin)
%
% arguments:  datarun - Standard datarun struct
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun with datarun.cones.mosaic struct set as in
%             make_mosaic_struct, preserving existing datarun.cones.mosaic
%             fields
%
%
% optional params, as in make_mosaic_struct
%
% 2011-05 phli, wrapper in standard datarun style and also to respect existing fields
%
%


mosaic = make_mosaic_struct(datarun.cones.centers);

if isfield(datarun.cones, 'mosaic')
    if isstruct(datarun.cones.mosaic)
        datarun.cones.mosaic = setstructfields(datarun.cones.mosaic, mosaic);
    else
        error('datarun.cones.mosaic exists but is not a struct');
    end
else
    datarun.cones.mosaic = mosaic;
end