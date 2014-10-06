function hexArray = elecHexArray(hexCoords)

if nargin < 1
    hexCoords = elecHexCoords();
end

qsize = max(hexCoords(:,1));
rsize = max(hexCoords(:,2));

hexArray = zeros(rsize, qsize);

for elec=1:512
    hexArray(hexCoords(elec, 2), hexCoords(elec, 1)) = elec;
end
end