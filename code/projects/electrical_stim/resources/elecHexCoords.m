%% script for generating axial hexagonal coordinates for electrodes
function hexCoords = elecHexCoords()
hexCoords = zeros(512, 2);
temp = load([matlab_code_path() 'private/freddy/512elecpositions.mat']); % Find a more general location for this or call a different text file.
positions = temp.positions;

startq = 9;
r = 1;
for ypos = fliplr(-390:120:450)
    q = startq;
    for xpos = -915:60:945
        elecID = ismember(positions, [xpos ypos], 'rows');
        hexCoords(elecID, :) = [q r];
        q = q + 1;
    end
    ypos2 = ypos - 60;
    r = r + 1;
    startq = startq - 1;
    q = startq;
    for xpos = -945:60:915
        elecID = ismember(positions, [xpos ypos2], 'rows');
        hexCoords(elecID, :) = [q r];
        q = q + 1;
    end
    r = r + 1;
end
end