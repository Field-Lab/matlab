function plotElectricRF(electrodes, sensitivities, markerScale, pElec)


[xCoords yCoords] = getElectrodeCoords61();

figure('position', [100 100 200 200])
axes('position', [0.2 0.2 0.6 0.6])
hold on
for i = 1:length(electrodes)
    markerSize = abs(sensitivities(i))*markerScale;
    if sensitivities(i) > 0
        plot(xCoords(electrodes(i)), yCoords(electrodes(i)), 'ro', 'markerSize', markerSize)
    else
        plot(xCoords(electrodes(i)), yCoords(electrodes(i)), 'bo', 'MarkerSize', markerSize)
    end
    text(xCoords(electrodes(i))+2, yCoords(electrodes(i)), num2str(sensitivities(i)))
end
plot(xCoords(pElec), yCoords(pElec), 'k*')
hold off
axis equal