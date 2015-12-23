hexCoords = elecHexCoords();
t = cputime;
hexNeighborsForElec(1:512, hexCoords);
totalTime = cputime-t;
disp(['Neighbor calculation (slow): ' num2str(totalTime)]);


t = cputime;
hexNeighborsFast(1:512, hexCoords);
totalTime = cputime-t;
disp(['Neighbor calculation (fast): ' num2str(totalTime)]);