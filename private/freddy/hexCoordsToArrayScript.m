function hexArray = elecHexArray()
hexCoords = elecHexCoords();

hexArray = zeros(16, 40);

for elec=1:512
    hexArray(hexCoords(elec, 2), hexCoords(elec, 1)) = elec;
end
end