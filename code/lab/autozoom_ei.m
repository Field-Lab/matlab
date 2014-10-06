function bounds = autozoom_ei(datarun, cellid, varargin)
% AUTOZOOM_EI 
% usage: bounds = autozoom_ei(datarun, cellid, varargin)
%
%

opts = inputParser();
opts.addParamValue('keep_aspect_ratio', true);
opts.addParamValue('axes', gca);
opts.addParamValue('aspect_ratio', 1);
opts.addParamValue('padding', []);
opts.addParamValue('padding_factor', 5);
opts.addParamValue('positions', datarun.ei.position);
opts.addParamValue('stack', []);
opts.parse(varargin{:});
opts = opts.Results;


if ~isempty(opts.stack)
    switch class(opts.stack)
        case 'cell'
            stack = datarun.stacks{opts.stack{:}};
        case 'struct'
            stack = opts.stack;
    end
    
    ai = parse_stack_index('array');
    array_tform = stack.tforms_inv{ai{:}};
    opts.positions = tformfwd(array_tform, opts.positions);
end


java_ei = datarun.ei.java_ei;
% Old style: maxelectrode = java_ei.getMaxElectrode(java_ei.getImage(cellid));
maxelectrode = java_ei.getMaxElectrode(cellid);
position = opts.positions(maxelectrode,:);


% If no specific padding given, calculate based on padding factor which is
% a multiple of the minimum electrode spacing
if isempty(opts.padding)
    spacing = ipdm(opts.positions, 'Subset', 'SmallestFew', 'Limit', 1);
    spacing = full(spacing(spacing > 0));
    opts.padding = opts.padding_factor * spacing;
end


if opts.keep_aspect_ratio
    curraxis = axis(opts.axes);
    xdiff = curraxis(2) - curraxis(1);
    ydiff = curraxis(4) - curraxis(3);
    opts.aspect_ratio = ydiff / xdiff;
end


xdiff = opts.padding;
ydiff = opts.padding * opts.aspect_ratio;
bounds = [position(1)-xdiff position(1)+xdiff position(2)-ydiff position(2)+ydiff];


if nargout < 1
    axis(opts.axes, bounds);
    clear bounds;
end
