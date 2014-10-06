komp='pokoj109';
if komp=='pokoj109'
    RawDataPath='D:\Home\Data\slices\2010-09-14-0\data002'; %pokoj 109
    paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\Home\Pawel\analysis\2010-09-14-0\data002min009new2\2010-09-14-0\data002\data002.params');
    neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\Home\Pawel\analysis\2010-09-14-0\data002min009new2\2010-09-14-0\data002\data002.neurons');
    dane_path='D:\Home\Pawel\analysis\2010-09-14-0\data002min009new2\dane';
else
    RawDataPath='C:\pawel\nauka\dane\2010-09-14-0\data002min009'; %laptop
    paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\pawel\nauka\analiza\slices\2010-09-14-0\data002min009new2\2010-09-14-0\data002\data002.params');
    neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\pawel\nauka\analiza\slices\2010-09-14-0\data002min009new2\2010-09-14-0\data002\data002.neurons');
    dane_path='C:\pawel\nauka\analiza\slices\2010-09-14-0\data002min009new2\dane;'
end
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(RawDataPath); 

idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=5241;
SpikeTimes = neuronFile.getSpikeTimes(NeuronID)';

cd(dane_path)%C:\pawel\nauka\analiza\slices\2010-09-14-0\data002min009new2\dane;
fid1=fopen(['ID=' num2str(NeuronID)],'r','ieee-le');
dane=fread(fid1,'int32');
fclose(fid1);
ElectrodesToCancel=dane(3:4:length(dane));
Latencies=dane(4:4:length(dane));

[PrimaryElectrode,Spikes,EI]=NS512_FindPrimaryElectrodeWithARtifact(RawDataPath,SpikeTimes,ElectrodesToCancel,Latencies);
PrimaryElectrode
figure(11)
plot(Spikes')