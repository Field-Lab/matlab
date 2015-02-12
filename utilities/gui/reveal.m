function reveal(h)
% 2011-07 phli

if nargin < 1
    h = gca;
end

f = getfig(h);
figure(f);
drawnow;