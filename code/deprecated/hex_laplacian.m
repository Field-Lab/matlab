function lv = hex_laplacian(v, neighbors)
% HEX_LAPLACIAN     Calculate the Laplacian on discretely sampled hexagonally spaced values
% usage: lv = hex_laplacian(v, neighbors)
%
% 2012-03 phli
%


lv = zeros(size(v));
for i = 1:length(neighbors)
    n = neighbors{i};
    lv(i) = length(n)*v(i) - sum(v(n));
end

lv = lv .* sqrt(1/3);