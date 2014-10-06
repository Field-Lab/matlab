function [timecourse, params] = time_course_from_sta(sta, sig_stixels, params)
% time_course_from_sta     get mean timecourse of significant stixels in a single STA
%
% usage:  [timecourse, tc_color] = time_course_from_sta(sta, sig_stixels, params)
%
% arguments:      sta - 4-d matrix of the STA (height, width, color, frames)
%         sig_stixels - 2D boolean matrix identifying significant stixels
%              params - struct of optional parameters (see below)
%
% outputs: timecourse - FxC matrix giving the timecourse
%                           F is the number of frames
%                           C is the number of color channels in the STA
%              extras - struct of more information (see below)
%              params - struct of params which were used
%
% optional fields in params, their default values, and what they specify:
%
%  color        []       If empty, returns the time course for all colors
%                        in STA. If it equals 1, 2, or 3 it returns only
%                        the red, green, or blue time course respectively
%
% gauthier, 2008-03
% edits:
%   gdf, 2008-09, now function must be handed sig_stixels
%

% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.color = [];
defaults.plot = false;
defaults.cell_name = [];

% combine user and default parameters
params = default_params( defaults, params);

% check inputs arguments
if nargin < 2
    error('check usage of function')
end
if ~islogical(sig_stixels)
    error('sig_stixels must be a logical array (boolean)')
end

% check defaults color is in range
if ~isempty(params.color) && params.color ~= 1 && params.color ~= 2 && params.color ~= 3
    error('specified color to return is not 1 (red), 2 (green), 3 (blue), or [] (all)')
end

% if sig_stixels is all 0s, return empty timecourse
if ~any(any(sig_stixels))
    timecourse = [];
    return
end


% COMPUTE TIME COURSE

% initialize timecourse variable
temp_timecourse = zeros(size(sta,4),size(sta,3));

% change to index values in the reshaped_sta (below)
sig_stixel_indices = find(reshape(sig_stixels,[],1));

% get mean timecourse of these stixels 
for cc = 1:size(sta,3)
% reshape the STA so that i = stixel identity, j = timecourse
    reshaped_sta = reshape(sta(:,:,cc,:),[],size(sta,4));
    temp_timecourse(:,cc) = squeeze(mean(reshaped_sta(sig_stixel_indices,:),1));
end

if params.color == 1    %red
    timecourse = temp_timecourse(:,1);
elseif params.color == 2    %green
    timecourse = temp_timecourse(:,2);
elseif params.color == 3    %blue
    timecourse = temp_timecourse(:,3);
else    %red, green, and blue, or black and white
    timecourse = temp_timecourse;
end


