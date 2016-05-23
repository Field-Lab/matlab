DataPath = 'Z:/Data/scan_new';
PatternNumber = 1; 
MovieNumber = 51;
Resolution = 10;
BadElectrodes = [353 354 420];

[line, v0, ve, vep, edges, func] = getPropagDirection (DataPath,PatternNumber,MovieNumber,BadElectrodes);

[Plots,Distances] = getCompleteIntersectionSet (DataPath,PatternNumber,MovieNumber,Resolution,BadElectrodes);

xx = -945:10:945;
yy = -450:10:450;
[qx,qy] = meshgrid(xx,yy);

z = func(qx,qy);
surf(qx,qy,z);

y = polyval(line,xx);
plot(xx,y);
axis([-945 945 -450 450])