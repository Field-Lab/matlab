function PrimaryNeurons=NS512_IdentifyPrimaryNeurons_v2(NeuronFilePath,DuplicatesFilePath)
%NeuronFilePath - like:
%'D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.neurons'
%DuplicatesFilePath - like:
%like:'D:\home\pawel\analysis\slices\2013\2013-12-12-3-PH\duplicates.txt'

neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(NeuronFilePath);
idList = neuronFile.getIDList();

fid=fopen(DuplicatesFilePath);
C=textscan(fid,'%d%s%s%d%s%f32');
D=C{1};
Duplicates=unique(D);
fclose(fid);

PrimaryNeurons=setdiff(idList,Duplicates);