function linkaxesandfit(axes, varargin)
% LINKAXESANDFIT    Link axes, keeping largest YLims
% usage: linkaxesandfit(axes, opts)
%
% Currently only handles YLim, but could be expanded for other properties
%
% 2012-03 phli
%

opts = inputParser();
opts.addParamValue('YLim', true);
opts.parse(varargin{:});
opts = opts.Results;


if islogical(opts.YLim) && opts.YLim
    ylims = zeros(numel(axes),2);
    for i = 1:numel(axes), ylims(i,:) = get(axes(i), 'YLim'); end
    opts.YLim = [min(ylims(:,1)) max(ylims(:,2))];
end


linkaxes(axes);
set(axes(1), 'YLim', opts.YLim);


for i = 1:numel(axes)
    set(axes(i), 'YTickMode', 'auto'); 
end