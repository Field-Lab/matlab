% generates a couple of figures that illustrate the identification of axonal signals for fitting a
% line to the axon (used for ARVO 2009 poster)

pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data006/data006.ei';

neuronID = 151;
[xCoords yCoords] = getElectrodeCoords61();

% get neuron's average spike waveforms from VISION analysis
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
ei = eiFile.getImage(neuronID); %gets ei data for neuron, storing the information as a 3D array: 2 x nElectrodes x nSamples
eiFile.close()

% calculate maximum waveform value on each electrode (absolute value)
eiAmps = zeros(64, 1);
for i = 1:64
    if ~(i==9||i==25||i==57)
        eiAmps(i) = max(max(abs(ei(1,i+1,:))));
    end
end

eiAmpsNormalized = eiAmps/max(eiAmps);


arrayWidth = max(xCoords)*2*1.2;
arrayHeight = max(yCoords)*2*1.2;

figure('position', [100 100 400 350])

for i = 1:64
    if ~(i==9||i==25||i==57)
        axes('position', [xCoords(i)/arrayWidth + 0.45, yCoords(i)/arrayHeight + 0.45, 0.1, 0.1])
        plot(squeeze(ei(1,i+1,:)), 'k')
        set(gca, 'ylim', [-10 8], 'xlim', [1 81], 'xtick', [], 'ytick', [])

        if any(i == [22 15 24 26 28 29 32 33 36 34 35 38 43 39 37 40 42])
            set(gca, 'box', 'on')
        else
            axis off
        end            
    end
end

figure('position', [100 100 400 350])
hold on
axis equal

for i = 1:64
    if ~(i==9||i==25||i==57)
        if eiAmpsNormalized(i)>.02
            if any(i == [22 15 24 26 28 29 32 33 36 34 35 38 43 39 37 40 42])
                plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*40), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
            else
                plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*40), 'MarkerFaceColor', 'k')
            end
        end
    end
end

hold off
