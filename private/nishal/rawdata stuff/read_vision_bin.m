
if 1% read single .bin data file 
    
    % open raw data file
    rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile('/Volumes/Data/2005-04-26-1/data006/data006000.bin');

    % gets data from sample 0 to 20000 (first second) across all electrodes
    samplingRate = 20000;
    data = rawFile.getData(0, samplingRate);

    % Stimulus TTLs (most times every 100 frames) on java electrode position 0 (1 in matlab)
    figure; 
    plot(data(:,1:3));

    rawFile.close();
end
  



if 0%alternate method reads across .bin data files
    
    % open raw data file 
    rawFile = edu.ucsc.neurobiology.vision.io.RawDataWrapper('/Volumes/Data/2005-04-26-1/data006');

    samplingRate = 20000;

    %set up arguments
    sample = 1000; %sample to start reading from
    nSamples = 1*samplingRate; % get one second, starting from sample
    nElectrodes = 513; % ALL electrodes in the dataset

    %get the raw data
    data = rawFile.getData(sample, nSamples, nElectrodes);

    % Stimulus TTLs (most times every 100 frames) on java electrode position 0 (1 in matlab)
    figure; 
    plot(data(:,1:4));

    rawFile.close();
end


%%

if 1% get electrode positions
    
    array_id=501; %placeholder for all rectangular arrays
    
    % get electrodeMap object
    electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(array_id);

    % get location of each electrode
    num_positions = electrodeMap.getNumberOfElectrodes - 1; % don't include TTLs
    positions = zeros(num_positions,2);
    for pp = 1:num_positions
        positions(pp,1) = electrodeMap.getXPosition(pp);
        positions(pp,2) = electrodeMap.getYPosition(pp);
    end
    
    %positions in micrometer
    plot(positions(:,1),positions(:,2),'.')    
end


hold on
list=[31,23,39,35,27,36,28];
plot(positions(list,1),positions(list,2),'o');



             
                