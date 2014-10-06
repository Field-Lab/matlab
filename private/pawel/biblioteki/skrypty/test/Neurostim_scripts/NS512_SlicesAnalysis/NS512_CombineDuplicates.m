neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.neurons');
DuplicatesFilePath='D:\home\pawel\analysis\slices\2013\2013-12-12-3-PH\duplicates.txt';

idList = neuronFile.getIDList();

a=importdata(DuplicatesFilePath);
sa=size(a.textdata);

VectorOfDuplicates=zeros(1,sa(1));
for i=1:sa(1) %unfortunately, the loop is required    
    VectorOfDuplicates(i)=str2num(a.textdata{i,1});
end

PrimaryNeurons=setdiff(idList,VectorOfDuplicates);

GoodPrimaryNeurons=PrimaryNeurons; % here, ultimately remove 'bad neurons' that should be saved in a separate text file,
%based on manual analysis in vision of the primary neurons!

PrimaryNeuronID=289;
Duplicates=[];
for i=1:sa(1) %unfortunately, the loop is required    
    PrimaryNeuron=str2num(a.textdata{i,4});
    if PrimaryNeuron==PrimaryNeuronID
        Duplicates=[Duplicates str2num(a.textdata{i,1})];
    end            
end