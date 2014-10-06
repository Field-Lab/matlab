ParamsFilePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.params';
NeuronFilePath = 'D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.neurons';
MovieFilePath='I:\analysis\slices\2013-12-12-3-PH\movie005';
DuplicatesFilePath='D:\home\pawel\analysis\slices\2013\2013-12-12-3-PH\duplicates.txt';

spikeTimes=NS512_SpikeTimesToStimulationParameters(PrimaryNeurons(1),NeuronFilePath,MovieFilePath,DuplicatesFilePath);