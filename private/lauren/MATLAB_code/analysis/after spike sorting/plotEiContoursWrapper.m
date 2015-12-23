pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data006/data006.ei';
neuronID = 151;
%pElec = 60;

% get neuron's average spike waveforms from VISION analysis
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
ei = eiFile.getImage(neuronID); %gets ei data for neuron, storing the information as a 3D array: 2 x nElectrodes x nSamples
clear eiFile

% calculate maximum waveform value on each electrode (absolute value)
eiAmps = zeros(64, 1);
for j = 1:64
    if ~(j==9||j==25||j==57)
        eiAmps(j) = max(max(abs(ei(1,j+1,:))));
    end
end

[xCoords yCoords] = getElectrodeCoords61();
eiAmps = eiAmps/max(eiAmps);
figure
hex_contour(xCoords, yCoords, eiAmps, 8, 'fig_or_axes', gca, 'contourSpacing', 'linear', 'plotCoords', false);
hold on

%primary electrode
%plot(xCoords(pElec), yCoords(pElec), 'k*')

%other electrodes
for i = 1:64
    %if ~any([9 25 57 pElec] == i)
        plot(xCoords(i), yCoords(i), 'ok','MarkerSize', 2, 'MarkerFaceColor', 'k')
    %end
end


%array outline
plot([0 8.6603 8.6603 0 -8.6603 -8.6603 0], [10 5 -5 -10 -5 5 10], 'k-')

set(gca, 'XLim', [-10 10], 'YLim', [-11 11])
axis equal
axis off
hold off