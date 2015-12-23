function check_alignment(stacks, input_index, base_index, varargin)
% CHECK_ALIGNMENT    Plot control points to confirm coordinate transforms
% usage: check_alignment(stacks, input_index, base_index, opts)
%
% opts: input_slice     1
%       base_slice      1
%
% 2010-08 phli
%

opts = inputParser;
opts.addParamValue('input_slice', 1);
opts.addParamValue('base_slice',  1);
opts.parse(varargin{:});
opts = opts.Results;

input_index = parse_stack_index(input_index);
base_index  = parse_stack_index(base_index);

istack = get_stack(stacks, input_index);
bstack = get_stack(stacks, base_index);

ips = get_input_points(istack, base_index);
bps = get_base_points( istack, base_index);


figure;
imshow(get_slice(istack, opts.input_slice));
hold on;
if ~isempty(ips), plot_xy(ips, '*'); end
if ~isempty(bps)
    tps = tformfwd(istack.tforms_inv{base_index{:}}, bps);
    plot_xy(tps, '*r');
    
    if ~isempty(ips)
        plot([ips(:,1) tps(:,1)]', [ips(:,2) tps(:,2)]', 'y');
    end
end


figure;
imshow(get_slice(bstack, opts.base_slice));
hold on;
if ~isempty(bps), plot_xy(bps, '*'); end
if ~isempty(ips)
    tps = tformfwd(istack.tforms{base_index{:}}, ips);
    plot_xy(tps, '*r');

    if ~isempty(bps)
        plot([bps(:,1) tps(:,1)]', [bps(:,2) tps(:,2)]', 'y');
    end
end
