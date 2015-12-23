% class: edu.ucsc.neurobiology.vision.io.NeuronFile
%    methods of interest:
%       getHeader()
%       getTTLTimes()
%       getSpikeTimes(neuronid)
      
% add the vision jar file to the matlab java path 
visionWritePath = '/Users/grosberg/code/write data GG/WriteDataFile.jar';
visionPath = '/Applications/Vision.app/Contents/Resources/Java/Vision.jar'; 
javaaddpath(visionWritePath);
javaaddpath(visionPath);

% for the neuron times, once you've got the neuron file, it's pretty simple
% add the vision jar file to the matlab java path 
pathToNeuronFile = '/Volumes/Analysis/2014-08-13-0/data003visiontest/2014-08-13-0/data012/data012.neurons';
neuronFile  = edu.ucsc.neurobiology.vision.io.NeuronFile(pathToNeuronFile);
header = neuronFile.getHeader();
spikes = neuronFile.getSpikeTimes(1396);

eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile('/Volumes/Analysis/2014-08-13-0/data003visiontest/2014-08-13-0/data012/data012.ei');
% will instantiate a neuron file object which you can use to interface to the file
neuronFile.getNeuronIDs() %will return a list of neuron IDs
neuronFile.getSpikeTimes(neuronid) %will return an array of spike times in samples
neuronFile.getTTLTimes() %will return the time of the triggers seen on electrode 0
neuronFile.getHeader().getSamplingFrequency() %will give you the sampling frequency for this data set which should usually be 20kHz


%% Try using other functions. 
piece = '/Volumes/Analysis/2014-08-13-0/data003visiontest/2014-08-13-0/'; 
datarun = load_data([piece 'data012/data012']);
datarun = load_neurons(datarun);