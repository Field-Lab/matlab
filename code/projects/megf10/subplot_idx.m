function [idx, xx, yy] = subplot_idx(x, y)

xx = 4*x-1; yy = 4*y-1;
idx = zeros(x*y, 9);
for i = 1:y
    idx(i, :) = [yy+2 yy+3 3 2 1 yy+1 2*yy+1 2*yy+2 2*yy+3] + 4*(i-1);
end

if x > 1
    for i = 2:x
        idx(y*(i-1)+1:y*i, :) = idx(1:y, :) + yy*4*(i-1);
    end
end

end

