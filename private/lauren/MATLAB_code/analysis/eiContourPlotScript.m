
%pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh/data001-lh.ei';
%neuronIDs = [677 676 754];

pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh/data001-lh.ei';
neuronIDs = [754];
%neuronIDs = [886 18 91 857 858 782 872 931];
%pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh/data001-lh.ei';

neuronIDs = 559;
pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data007-NW/data007-NW.ei';




% figure
% axes
% plotEi61(pathToEi, neuronIDs)

[xCoords yCoords] = getElectrodeCoords61();

CsOM = cell(length(neuronIDs), 1);
for i = 1:length(neuronIDs)
    neuronID = neuronIDs(i);
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
    eiAmps = eiAmps/max(eiAmps);

%     CsOM{i} = getCOM(eiAmps, 0.2, 1);
    
    contours = hex_contour(xCoords, yCoords, eiAmps, 8, 'fig_or_axes', 0, 'contourSpacing', 'linear');
    title(['neuron' num2str(neuronIDs(i))])
end

%array outline
hold on
plot([0 8.6603 8.6603 0 -8.6603 -8.6603 0], [10 5 -5 -10 -5 5 10], 'k-')

set(gca, 'XLim', [-10 10], 'YLim', [-11 11])
axis equal
axis off

% %plot of electrode positions
% figure
% hold on
% for i = 1:64
%     plot(xCoords(i), yCoords(i), 'ok','MarkerSize', 2, 'MarkerFaceColor', 'k')
% end
% plot([0 8.6603 8.6603 0 -8.6603 -8.6603 0], [10 5 -5 -10 -5 5 10], 'k-')
% hold off
% 
% set(gca, 'XLim', [-10 10], 'YLim', [-11 11])
% axis equal
% axis off
