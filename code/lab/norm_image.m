function im = norm_image(sta_frame, params)
% norm_image     adjust contrast to plot STA frame.  maps 0 to 0.5, and scales so all values are between 0 and 1.
%
% usage:  im = norm_image(sta_frame)
%
% arguments:  sta_frame - 3-d matrix (x,y,color)
%                params - struct of optional parameters (see below)
%
% outputs:           im - image for plotting in Matlab format
%                           all values between 0 and 1
%
%
% optional fields in params, their default values, and what they specify:
%
% max               []  	maximum contrast.  if empty, determined automatically
% rgb               true    if true, will expand a one-color summary frame to have 
%                               3 (identical) color channels
%
%


% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.max = [];
defaults.rgb = true;

% combine user and default parameters
params = default_params( defaults, params);


% BODY OF THE FUNCTION

% if sta_frame is empty or a scalar, return empty matrix
if isempty(sta_frame) || numel(sta_frame) == 1
    im = [];
    return
end


% convert to double
sta_frame = double(sta_frame);

% get maximum contrast
if isempty(params.max)
    params.max = max(max(max(max(abs(sta_frame)))));
end

% normalize so that maximum contrast = 1
im = sta_frame/params.max;

% clip values outside this range
im(im < -1) = -1;
im(im >  1) =  1;

% normalize so that -1 -> 0, 0 -> 0.5, 1 -> 1
im = 0.5 + im/2;

% ensure there are three color dimensions
if params.rgb && (size(im,3) == 1)
    im = repmat(im,[1 1 3]);
end
