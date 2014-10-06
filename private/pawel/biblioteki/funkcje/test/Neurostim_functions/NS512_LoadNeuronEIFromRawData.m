function NeuronEI=NS512_LoadNeuronEIFromRawData(RawDataPath,paramsFilePath,neuronFilePath,NeuronID,MaxSpikeNumber,Length,Offset);
%Ta funkcja laduje EI konkretnego neuronu, uzywajac funkcji
%NS512_LoadEIFromRawData. Argumenty wejsciowe:
%MaxSpikeNumber - jesli dla danego neuronu jest np. 20 tys. spikow, mozemy
%nie chciec ladowac wszystkich dla danego EI. Ten parametr okresla, ile
%spikw powinno byc wzietych. Jesli np. wartosc tego parametru to 100,
%wowczas konkretne numery spikow sa potem brane do EI to: 1, 201, 401, ...

paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile(paramsFilePath);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronFilePath);

spikeTimesAll= neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

SpikeTimesStep=floor(length(spikeTimesAll)/MaxSpikeNumber);

spikeTimes=spikeTimesAll(1:SpikeTimesStep:1+SpikeTimesStep*(MaxSpikeNumber-1));

NeuronEI=NS512_LoadEIFromRawData(RawDataPath,spikeTimes,Length,Offset);