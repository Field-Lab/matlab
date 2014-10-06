function PrimaryNeurons=NS512_IdentifyPrimaryNeurons(NeuronFilePath,DuplicatesFilePath)
%NeuronFilePath - like:
%'D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.neurons'
%DuplicatesFilePath - like:
%like:'D:\home\pawel\analysis\slices\2013\2013-12-12-3-PH\duplicates.txt'

neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(NeuronFilePath);
idList = neuronFile.getIDList();

a=importdata(DuplicatesFilePath);
sa=size(a.textdata)

VectorOfDuplicates=zeros(1,sa(1));
for i=2:sa(1) %unfortunately, the loop is required    
    a.textdata{i,1}
    VectorOfDuplicates(i)=str2num(a.textdata{i,1});
end

PrimaryNeurons=setdiff(idList,VectorOfDuplicates);