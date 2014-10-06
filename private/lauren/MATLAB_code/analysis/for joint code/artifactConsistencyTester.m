cd /netapp/snle/home/lhruby/Desktop/testFigs/artifactDifferences/2008-08-26-0-data006

dataPath1 = '/Volumes/Lee/Analysis/Lauren/';
dataPath2_1 = '2008-08-26-0';
dataPath2_2 = 'data006';
dataPath2 = [dataPath2_1, '/', dataPath2_2];
dataPath = [dataPath1 dataPath2];

artifactPath1 = '/Volumes/Lee/Analysis/Lauren/';
artifactPath2 = '2008-08-26-0/data011';
artifactPath = [artifactPath1 artifactPath2];

%clusterPath = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/clusters011artifacts';
clusterPath = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/clusters006artifacts';

%dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data001/';
%artifactPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data013/';

%dataPath = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data011';
%artifactPath = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data023';

centerChannel = 21;
patternNumber = 21;
PW = 100;
doPlot = 0;

patternNumbers = patternNumber;
%movieNumbers = 1:26;
movieNumbers = 1:1:66;

artPatternNumbers = patternNumber;
%artMovieNumbers = 2:3:74;
artMovieNumbers = 1:1:66;


nMovies = length(movieNumbers);
nPatterns = length(patternNumbers);


%% making plots

lineColors = hsv(26);

figure(1)
clf
figure(2)
clf
figure(3)
clf

amps = zeros(1, nMovies);
meanArtifact = cell(1, nMovies);
meanDataTrace = cell(1, nMovies);
meanDifference = cell(1, nMovies);

for i = 1:nMovies
    clusterIndex = NS_ReadClusterFile(clusterPath, movieNumbers(i), patternNumbers(1));
    
    amps(i) = getStimAmps(dataPath, patternNumbers(1), movieNumbers(i));
    
    dataTraces=NS_ReadPreprocessedData(dataPath, '', 0, patternNumbers(1), movieNumbers(i), 99999);
    artifactTraces = NS_ReadPreprocessedData(artifactPath, '', 0, artPatternNumbers(1), artMovieNumbers(i), 99999);
    
    meanArtifact{i} = mean(artifactTraces(:,centerChannel,:), 1);
    meanDataTrace{i} = mean(dataTraces(clusterIndex==1, centerChannel, :),1);
    meanDifference{i} = squeeze(meanDataTrace{i} - meanArtifact{i});
    
    if doPlot
        figure(2)
        hold on
        current = plot(squeeze(meanDataTrace{i}));
        set(findobj(current,'Type','line'),'Color',lineColors(i,:))
        hold off

        figure(1)
        hold on
        current = plot(squeeze(meanArtifact{i}));
        set(findobj(current,'Type','line'),'Color',lineColors(i,:))
        hold off

        figure(3)
        hold on
        current = plot(meanDifference{i});
        set(findobj(current,'Type','line'),'Color',lineColors(i,:))
        hold off
    end
end

if doPlot
    hold off
    title(['mean difference between normal recording and TTX recording at different current amplitudes', 10,...
        'current range ', num2str(amps(1)), ' - ', num2str(amps(nMovies)), '\muA, PW: ', num2str(PW),  '\mus, electrode' num2str(centerChannel), 10,...
        'data:', dataPath2, ', pattern ', num2str(patternNumbers), 10,...
        'artifact:', artifactPath2, ', pattern', num2str(artPatternNumbers)])

    xlabel('samples')
    ylabel('DAQ units')
end
    
details.dataset = dataPath2;
details.artifactDataset = artifactPath2;
details.amps = amps;
details.patternNumber = patternNumbers(1);
details.dataMovieNumbers = movieNumbers;
details.artifactMovieNumbers = artMovieNumbers(1:length(movieNumbers));

details.description = ['Data traces with clear spikes (up to sample 35) were removed before taking their average.', 10,...
    'Included movie numbers represent stimulation amplitudes that elicited no obvious spikes in at least a subset',...
    ' of the data traces.'];


save([dataPath2_1 '_' dataPath2_2 '_p' num2str(patternNumbers(1))], 'meanArtifact', 'meanDataTrace', 'meanDifference', 'details')


