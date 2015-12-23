function lv = lazyhex_laplacian(v, neighbor_struct)
% LAZYHEX_LAPLACIAN     Calculate the Laplacian on discretely sampled lazy-hexagonally spaced values
% usage: lv = hex_laplacian(v, neighbor_struct)
%
% Our arrays are lazy-hexagonal; pitch L along a row, pitch L between rows,
% offset L/2 for adjacent rows.
%
% 2012-03 phli
%


lv = zeros(size(v));
for i = 1:length(neighbor_struct)
    n = neighbor_struct(i);
    lv(i) =         2*length(n.acrosslines)*v(i) - 2*sum(v(n.acrosslines));
    lv(i) = lv(i) + 3*length(n.inline)     *v(i) - 3*sum(v(n.inline));
end

lv = lv ./ 14;