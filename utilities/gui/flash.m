function flash(h)
% 2011-07 phli

if nargin < 1
    h = gca;
end

oldvis = get(h, 'visible');

for i = 1:10
    proptoggle(h, 'visible');
    drawnow;
    pause(0.07);
end

if iscell(oldvis)
    for i = 1:length(h)
        val = oldvis(i);
        if iscell(val), val = val{1}; end
        set(h(i), 'visible', val);
    end
else
    set(h, 'visible', oldvis);
end