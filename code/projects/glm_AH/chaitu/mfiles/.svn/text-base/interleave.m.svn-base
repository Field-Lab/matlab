% Interleave 2 vectors

function c = interleave(a,b)


c = zeros(length(a)+length(b),1);

if (length(a) <= length(b))
    vmin = a;
    vmax = b;
else
    vmin = b;
    vmax = a;
end
1;

c(1:2:(2*length(vmin)-1)) = vmin;
c(2:2:(2*length(vmin))) = vmax(1:length(vmin));
c(2*length(vmin)+1:end) = vmax(length(vmin)+1:end);
