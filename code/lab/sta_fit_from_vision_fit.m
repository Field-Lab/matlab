function sta_fit = sta_fit_from_vision_fit(datarun,vision_fit)
% sta_fit_from_vision_fit     convert vision fit to sta fit in standard form
%
% usage:  sta_fit = sta_fit_from_vision_fit(datarun,vision_fit)
%
% arguments:     datarun - datarun struct
%             vision_fit - struct with fields mean, sd, angle
%
% outputs:       sta_fit - sta fit in standard form (mean, sd, angle, surround_sd_scale, surround_scale)
%
%
% 2009-04  gauthier
%
% NOTE: currently only spatial fit is supported, NOT time or color
%


% error check
if ~isfield(datarun,'stimulus') || ~isfield(datarun.stimulus,'field_height')
    error('cannot load vision fit because stimulus parameters are not known')
end


% BODY OF THE FUNCTION

% initalize output
sta_fit = struct;

% if input is empty, return empty
if isempty(fieldnames(vision_fit));return;end

% convert center point
sta_fit.mean = [vision_fit.mean(1) (datarun.stimulus.field_height - vision_fit.mean(2))] + 0.5;

% copy sd and angle
sta_fit.sd = vision_fit.sd;
sta_fit.angle = vision_fit.angle;

% fill in surround with empty values
sta_fit.surround_sd_scale = 1;
sta_fit.surround_scale = 0;