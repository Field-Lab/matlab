function plot_electrodes_over_slice(datarun, stack_spec, varargin)
% PLOT_ELECTRODES_OVER_SLICE    Plot the electrode locations over the stack slice image
% usage: plot_electrodes_over_slice(datarun, stack_spec, opts)
%
% opts: slice   1       Which slice to plot over
%       text    true    Whether to plot text electrode numbers
%       color   'b'     Color to use for plotting
%       plot    'o'     Plot command to use for plotting electrodes.  If
%                       set to empty, electrodes will not plot (i.e. you
%                       can have just the electrode numbers text)
%
% 2010-09 phli
%

opts = inputParser;
opts.addParamValue('slice', 1);
opts.addParamValue('text', true);
opts.addParamValue('plot', 'o');
opts.addParamValue('color', 'b');
opts.parse(varargin{:});
opts = opts.Results;

stack_index = parse_stack_index(stack_spec);
stack = get_stack(datarun, stack_index);
slice = get_slice(stack, opts.slice);

positions = datarun.ei.position;
positions = tformfwd(stack.tforms_inv{1}, positions);

imshow(slice);
hold on;

for i = 1:length(positions)
    if opts.plot
        plot(positions(i,1), positions(i,2), opts.plot, 'Color', opts.color);
    end
    
    if opts.text
        text(positions(i,1), positions(i,2), num2str(i), 'Color', opts.color);
    end
end