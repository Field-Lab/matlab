oldneuronfilepath = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data000/data000.neurons';
newneuronfilepath = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data000-debug/data000-debug.neurons';
rawdatafilepath = '/Volumes/Data/2013-08-19-3/data000';

neuronsPerElectrode= 10;

oldNeuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(oldneuronfilepath);
newNeuronFile = edu.ucsc.neurobiology.vision.matlab.ReadNeurons(newneuronfilepath, oldNeuronFile);

st = oldNeuronFile.getSpikeTimes(35);
for kk = 1:(length(st))
    newNeuronFile.addNeuron(floor(kk/neuronsPerElectrode)+1, ...
        mod(kk, neuronsPerElectrode), st(kk), 1);
end
% newNeuronFile.addNeuron(1, 0, st, length(st));

oldNeuronFile.close();
newNeuronFile.closeNeurons();