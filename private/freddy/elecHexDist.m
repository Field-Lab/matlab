function dist = elecHexDist(q1, r1, q2, r2)
    x1 = q1;
    z1 = r1;
    x2 = q2;
    z2 = r2;
    y1 = -(x1 + z1);
    y2 = -(x2 + z2);
    dist = (abs(x1 - x2) + abs(y1 - y2) + abs(z1 - z2)) / 2;
end