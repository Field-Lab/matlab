function [h, pos] = imarrow(varargin)
% IMARROW   Wrap for imline that adds some crappy directional indicators
%
% 2012-08 phli
%

h = imline(varargin{:});
c = get(h, 'Children');

starth = c(2);
endh   = c(1);
set(starth, 'Marker', 'o', 'MarkerSize', 7);
set(endh, 'Marker', 'p');

if nargout < 1, clear h; end
if nargout > 1
    toplineh = c(3);
    pos = [get(toplineh, 'XData'); get(toplineh, 'YData')]';
end