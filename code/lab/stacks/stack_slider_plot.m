function hout = stack_slider_plot(stack, handle, varargin)
% STACK_SLIDER_PLOT    Show an image stack with a GUI slider to look through sections
% usage: handle = stack_slider_plot(stack, [handle])
%
% Stack can be either a 3D monochrome XYZ matrix, a 4D X-Y-RGB-Z, or a stack struct as described in proposal.rtf.
%
% Recommend usage without passing in handle.
%
% Essentially just a wrapper on GEN_SLIDER_PLOT.
%
% 2010-04 phli
%


if nargin < 2 || isempty(handle)
    handle = figure();
end

opts = inputParser();
opts.addParamValue('channels', []);
opts.parse(varargin{:});
opts = opts.Results;


% Convert mono 3D stack to 4D RGB if necessary (stack can be either raw data, or a struct as described in proposal.rtf)
if ~isstruct(stack) && ndims(stack) == 3
    stack = reshape(stack, [size(stack, 1) size(stack, 2) 1 size(stack, 3)]);
end

% Customize api functions depending on whether stack is raw data or a struct as described in proposal.rtf
if isstruct(stack)
    stack_length_func = @stack_length;
    get_section_func  = @(stck, i) pick_channels(load_section(stck, i, handle), opts.channels);
else
    if isempty(opts.channels), opts.channels = 1:size(stack,3); end
    stack_length_func = @(stck) size(stck, 4);
    get_section_func  = @(stck, i) pick_channels(stack(:,:,:,i), opts.channels);
end

fig = gen_slider_plot('data', stack, 'handle', handle, 'data_length_func', stack_length_func, 'get_data_func', get_section_func, 'plot_func', @plot_func);


if nargout > 0
    hout = fig;
end



function h = plot_func(imag, sval, data, f)
v = axis;

ax = getappdata(f, 'slider_axes');
cla(ax);

h = imshow(imag);
title(sprintf('Section %i', sval));
axis fill;

if any(v ~= [0 1 0 1])
    axis(v);
end


% Method for getting sections and caching them for struct stacks
function section = load_section(stack, isect, f)
stack = load_slices(stack, isect);
setappdata(f, 'data', stack);
section = get_slice(stack, isect);


function im = pick_channels(im, channels);
if isempty(channels), return; end

% if length(channels) == 1
%     im = im(:,:,channels);
%     return;
% end

num_orig_chan = size(im,3);
kill_chan = setdiff(1:num_orig_chan, channels);
im(:,:,kill_chan) = 0;