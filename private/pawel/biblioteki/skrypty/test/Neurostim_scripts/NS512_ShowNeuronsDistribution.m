electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
RawDataPath='J:\2011-07-01-1\data000'; %define path to raw data file
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\Home\Pawel\analysis\slices\2011-06-29-1\data000\2011-06-29-1\data000\data000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\Home\Pawel\analysis\slices\2011-06-29-1\data000\2011-06-29-1\data000\data000.neurons');

%rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 

idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';

figure(3)
indeks=1;
X=[];
Y=[];
for neuron=1:length(idList)
    NeuronID=idList(neuron)
    %SpikeTimes = neuronFile.getSpikeTimes(NeuronID)';    
    SeedEl = neuronFile.getNeuronIDElectrode(NeuronID);
    if length(SpikeTimes)>200
        %[PrimaryElectrode,Spikes,EI]=NS512_FindPrimaryElectrodeWithARtifact(RawDataPath,SpikeTimes,[],[]);        
        %SeedEl = neuronFile.getNeuronIDElectrode(Neurons(i));
        X(indeks)=electrodeMap.getXPosition(SeedEl);
        Y(indeks)=electrodeMap.getYPosition(SeedEl);    
        indeks=indeks+1;
    end
end

length(unique(X+Y*10000))
plot(X,Y,'bd');
axis([-1000 1000 -500 500]);