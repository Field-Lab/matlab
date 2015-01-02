function setaxesy(ax, varargin)
% SETAXESY   Set axes YLims large enough to fit all data
% usage: setaxesy(axes, opts)
%
% 2012-03 phli
%

opts = inputParser();
opts.addParamValue('YLim', true);
opts.addParamValue('tight', false);
opts.parse(varargin{:});
opts = opts.Results;


if islogical(opts.YLim) && opts.YLim
    ylims = zeros(numel(ax),2);
    for i = 1:numel(ax)
        if opts.tight, axis(ax(i), 'tight'); end
        ylims(i,:) = get(ax(i), 'YLim'); 
    end
    opts.YLim = [min(ylims(:,1)) max(ylims(:,2))];
end


for i = 1:numel(ax)
    set(ax(i), 'YLim', opts.YLim);
    set(ax(i), 'YTickMode', 'auto'); 
end