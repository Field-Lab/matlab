[xCoords yCoords] = getElectrodeCoords61();

pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
cellID = 227;


% get neuron's average spike waveforms from VISION analysis
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);

ei = eiFile.getImage(cellID); %gets ei data for neuron,  as a 3D array: 2 x nElectrodes x nSamples
eiAmps = zeros(64, 1);
for j = 1:64
    if ~(j==9||j==25||j==57)
        eiAmps(j) = max(max(abs(ei(1,j+1,:))));
    end
end
eiAmpsNormalized = eiAmps/max(eiAmps);

clear eiFile

figure
hold on
for i = 1:64
    if ~(i==4||i==9||i==25||i==57) %% dead or nonexistent electrodes
        if ~(i==16||i==14||i==11||i==10||i==13||i==18||i==12||i==8||i==19) % electrodes with non-axonal signal
            if eiAmpsNormalized(i)>.01
                plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*40), 'MarkerFaceColor', 'k')
            end
        %plot(squeeze(ei(1,i+1,:)))
        end
    end
end
hold off

set(gca, 'XLim', [-10.4633 10.4633], 'YLim', [-8 8])
axis off
hold off



%% 243
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
cellID = 243;
% get neuron's average spike waveforms from VISION analysis
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);

ei = eiFile.getImage(cellID); %gets ei data for neuron,  as a 3D array: 2 x nElectrodes x nSamples
eiAmps = zeros(64, 1);
for j = 1:64
    if ~(j==9||j==25||j==57)
        eiAmps(j) = max(max(abs(ei(1,j+1,:))));
    end
end
eiAmpsNormalized = eiAmps/max(eiAmps);

clear eiFile

figure
hold on
for i = 1:64
    if ~(i==4||i==9||i==25||i==57) %% dead or nonexistent electrodes
        if ~(i==14||i==17||i==12||i==15||i==22||i==24||i==21||i==16||i==13||i==41||i==19||i==18) % electrodes with non-axonal signal
            if eiAmpsNormalized(i)>.02
                plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*100), 'MarkerFaceColor', 'k')
            end
        %plot(squeeze(ei(1,i+1,:)))
        end
    end
end
hold off

set(gca, 'XLim', [-10.4633 10.4633], 'YLim', [-8 8])
axis off
hold off

%%
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data012/data012.ei';

cellID = 334;
% get neuron's average spike waveforms from VISION analysis
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);

ei = eiFile.getImage(cellID); %gets ei data for neuron,  as a 3D array: 2 x nElectrodes x nSamples
eiAmps = zeros(64, 1);
for j = 1:64
    if ~(j==9||j==25||j==57)
        eiAmps(j) = max(max(abs(ei(1,j+1,:))));
    end
end
eiAmpsNormalized = eiAmps/max(eiAmps);

clear eiFile

figure
hold on
for i = 1:64
    if ~(i==4||i==9||i==25||i==57) %% dead or nonexistent electrodes
        if ~(i==27||i==24||i==23||i==29||i==30||i==19||i==20||i==28||i==18||i==40||i==37||i==39||i==42||i==16) % electrodes with non-axonal signal
            if eiAmpsNormalized(i)>.01
                plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*50), 'MarkerFaceColor', 'k')
            end
        %plot(squeeze(ei(1,i+1,:)))
        end
    end
end
hold off

set(gca, 'XLim', [-10.4633 10.4633], 'YLim', [-8 8])
axis off
hold off


%% 558

pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data009/data009.ei';
cellID = 558;

% get neuron's average spike waveforms from VISION analysis
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);

ei = eiFile.getImage(cellID); %gets ei data for neuron,  as a 3D array: 2 x nElectrodes x nSamples
eiAmps = zeros(64, 1);
for j = 1:64
    if ~(j==9||j==25||j==57)
        eiAmps(j) = max(max(abs(ei(1,j+1,:))));
    end
end
eiAmpsNormalized = eiAmps/max(eiAmps);

clear eiFile

figure
hold on
for i = 1:64
    if ~(i==4||i==9||i==25||i==57) %% dead or nonexistent electrodes
        if ~(i==44||i==38||i==43||i==46||i==39||i==36||i==41||i==47||i==54||i==53||i==56||i==49) % electrodes with non-axonal signal
            if eiAmpsNormalized(i)>.01
                plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*100), 'MarkerFaceColor', 'k')
            end
        %plot(squeeze(ei(1,i+1,:)))
        end
    end
end
hold off

set(gca, 'XLim', [-10.4633 10.4633], 'YLim', [-8 8])
axis off
hold off

%% 886
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data012/data012.ei';
cellID = 886;

% get neuron's average spike waveforms from VISION analysis
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);

ei = eiFile.getImage(cellID); %gets ei data for neuron,  as a 3D array: 2 x nElectrodes x nSamples
eiAmps = zeros(64, 1);
for j = 1:64
    if ~(j==9||j==25||j==57)
        eiAmps(j) = max(max(abs(ei(1,j+1,:))));
    end
end
eiAmpsNormalized = eiAmps/max(eiAmps);

clear eiFile

figure
hold on
for i = 1:64
    if ~(i==4||i==9||i==25||i==57) %% dead or nonexistent electrodes
        if 1 % electrodes with non-axonal signal
            if eiAmpsNormalized(i)>.01
                plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*40), 'MarkerFaceColor', 'k')
            end
        %plot(squeeze(ei(1,i+1,:)))
        end
    end
end
hold off

set(gca, 'XLim', [-10.4633 10.4633], 'YLim', [-8 8])
axis off
hold off

