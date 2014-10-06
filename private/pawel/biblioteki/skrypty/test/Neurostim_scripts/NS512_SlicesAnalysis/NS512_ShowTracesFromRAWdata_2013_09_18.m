clear
NS_GlobalConstants=NS_GenerateGlobalConstants(512);

%1) Reading raw data
full_path='D:\Home\Data\slices\2010-09-14-0\data002'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 

paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\Vision_output\2010-09-14-0\data002\data002.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\Vision_output\2010-09-14-0\data002\data002.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
NeuronID=6653
patternID=13;
WhichMovie=12;
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
NS512_ShowTracesFromRAWdata_ciag_dalszy_2013_09_18;