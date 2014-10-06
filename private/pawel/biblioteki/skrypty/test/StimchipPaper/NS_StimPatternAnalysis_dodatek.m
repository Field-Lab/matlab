full_path='H:\2010-09-21-0\data002\'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 

%2) Reading some neuron information
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('E:\pawel\analysis\retina\2010-09-21-0\data002_old\data002000\data002000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('E:\pawel\analysis\retina\2010-09-21-0\data002_old\data002000\data002000.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
%NeuronsIDs=[227 391 406 541 616 736 901 856];

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
%StimElectrodes=[16 27 28 37 45 51 61 60];

Channels=[1:64];
Channels=[1:8 10:24 26:56 58:64];
N=500;
L=80;
spikes=zeros(length(Channels),L)

for i=1:length(NeuronsIDs)
    NeuronID=NeuronsIDs(i);
    spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
    
    for j=1:N
        t=spikeTimes(j+100);
        d0=rawFile.getData(t-20,L)';
        d1=d0(Channels+1,:); % add 1, since the first channel in the raw data is the TTL channel and not electrode 1
        spikes(:,:)=spikes(:,:)+double(d1);
        %plot(d1(61,:));
        %hold on;
    end
    spikes=spikes/N;
    %subplot(3,3,i);
    subplot('position',[0.003,0.956-i*0.125,0.08,0.11]);
    'sfg'
    for j=1:length(Channels)
        spikes(j,:)=spikes(j,:)-mean(spikes(j,1:2));
    end
    % Normalization:
    spikes=spikes/max(max(abs(spikes)))*65;
    h=NS512_ShowEIAsCircles(spikes,1,Channels,StimElectrodes(i),[-305 305],[-305 305]);    
    h=NS_AddFrameForArrayLayout(1,LineWidth);
    h=gca;
    set(h,'Visible','off');
end