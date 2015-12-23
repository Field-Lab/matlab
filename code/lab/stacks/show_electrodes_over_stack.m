function varargout = show_electrodes_over_stack(datarun, stack_index, varargin)
% SHOW_ELECTRODES_OVER_STACK    Plot electrode positions over stack image
% usage: show_electrodes_over_stack(datarun, stack_index)
%
% 2010-12 phli
%

opts = inputParser;
opts.addParamValue('slice_index', 1);
opts.addParamValue('include_outside', false);
opts.addParamValue('color', 'r');
opts.addParamValue('label', true);
opts.addParamValue('scale', 0.1);
opts.addParamValue('channels', []);
opts.parse(varargin{:});
opts = opts.Results;


stack_index = parse_stack_index(stack_index);
stack = get_stack(datarun, stack_index);
im = get_slice(stack, opts.slice_index);
if size(im,3) == 4, im = im(:,:,1:3); end
if ~isempty(opts.channels)
    chanim = zeros(size(im), class(im));
    chanim(:,:,opts.channels) = im(:,:,opts.channels);
    im = chanim;
end
imshow(im);
hold on;

positions = datarun.ei.position;
ai = parse_stack_index('array');
array_tform_inv = stack.tforms_inv{ai{:}};
positions = tformfwd(array_tform_inv, positions);

if ~opts.include_outside
    elecs_in_stack = electrodes_in_stack(datarun, stack_index);
    positions = positions(elecs_in_stack, :);
    labels = collect(num2cell(elecs_in_stack), @num2str);
else
    labels = {};
end

plot_electrodes(positions, 'scale', opts.scale, 'label', opts.label, 'labels', labels, 'pos_color', opts.color);
axis image