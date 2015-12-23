function id = electrodeIDForXY(positions, x, y)
id = 0;
for i = 1:size(positions, 1)
    if positions(i, 1) == x && positions(i, 2) == y
        id = i;
    end
end
end