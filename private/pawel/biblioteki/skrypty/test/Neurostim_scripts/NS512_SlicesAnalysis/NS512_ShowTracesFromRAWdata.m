clear
NS_GlobalConstants=NS_GenerateGlobalConstants(512);

%1) Reading raw data
full_path='C:\pawel\nauka\dane\2010-09-14-0\data002min009'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 

paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\pawel\nauka\analiza\slices\2010-09-14-0\data002minus009paramSninya\2010-09-14-0\data002min009\data002min009.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\pawel\nauka\analiza\slices\2010-09-14-0\data002minus009paramSninya\2010-09-14-0\data002min009\data002min009.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=1217
patternID=13;
WhichMovie=6;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

%Find EI
N=150;
L=60;
offset=50;
LBack=30;
spikes2=zeros(N,512,L); %spike,channel,sample
for i=1:N
    t=spikeTimes(offset+i);
    d0=rawFile.getData(t-LBack,L)';
    d1=d0([1:512]+1,:);
    spikes2(i,:,:)=d1;
end
EI=mean(spikes2);
Amplitudes=max(abs(EI),[],3); % max amplitude (any polarity) for each electrode for this EI
figure(3)
plot(Amplitudes)
MaxElectrodes=find(Amplitudes==max(Amplitudes));
PrimaryElectrode=MaxElectrodes(1);
gh=spikes2(:,PrimaryElectrode,:);
gh1=reshape(gh,N,L);
figure(1)
plot(gh1');
grid on
break
figure(2)
NS512_ShowTracesFromRAWdata_ciag_dalszy;