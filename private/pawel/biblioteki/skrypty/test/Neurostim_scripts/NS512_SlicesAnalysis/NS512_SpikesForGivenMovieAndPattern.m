NS_GlobalConstants=NS_GenerateGlobalConstants(512);
MovieData=NS_MovieData_GlobalPath('D:\Home\Data\slices\2010-09-14-0\movie002',100,NS_GlobalConstants); % of course define your own path and file name
FirstRepetitionTime=MovieData(3) % this will be always an integer multiplication of 10000
NumberOfRepetitions=MovieData(4) % repetition period is of course 10000

ElectrodeToRead=250;
NeuronID=6784;
MovieNumber=75;

TimeBrackets=[FirstRepetitionTime,FirstRepetitionTime+NumberOfRepetitions*10000]

RawDataPath='D:\Home\Data\slices\2010-09-14-0\data002'; %pokoj 109
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_June_2011\data002minus009paramSninya\2010-09-14-0\data002min009\data002min009.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_June_2011\data002minus009paramSninya\2010-09-14-0\data002min009\data002min009.neurons');
    
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';

%NS512_SpikeTimesToStimulationPatterns('D:\Home\Data\slices\2010-09-14-0\movie002',SpikeTime,MoviesBegins,NS_GlobalConstants);