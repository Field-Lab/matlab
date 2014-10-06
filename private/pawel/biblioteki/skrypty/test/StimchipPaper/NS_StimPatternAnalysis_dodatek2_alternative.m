full_path='F:\2010-09-21-0\data002\'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 

StimElectrodes=[16 27 37 45 54  51 60];
RecElectrodes=[16 27 37 45 47 51 58];
NeuronsIDs=[227 391 541 616 691 736 856];
NeuronsIDs=[227 391 541 616 691 736 856];
NeuronsIDs=[541 616 691 736 856 901];
%NeuronsIDs=[76 227 256 271 391 406];


%2) Reading some neuron information
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\pawel\nauka\analiza\data002000\data002000.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\pawel\nauka\analiza\data002000\\data002000.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
%NeuronsIDs=[227 391 406 541 616 736 901 856];

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
%StimElectrodes=[16 27 28 37 45 51 61 60];

Channels=[1:64];
Channels=[1:8 10:24 26:56 58:64];
N=500;
L=80;
spikes=zeros(length(Channels),L)

for i=5
    %length(NeuronsIDs)
    NeuronID=NeuronsIDs(ElectrodeOrder(i));
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
    subplot('position',[0.003,0.95-i*0.145,0.09,0.12]);
    'sfg'
    for j=1:length(Channels)
        spikes(j,:)=spikes(j,:)-mean(spikes(j,1:2));
    end
    % Normalization:
    spikes=spikes/max(max(abs(spikes)))*70;
    h=NS512_ShowEIAsCircles(spikes,1,Channels,StimElectrodes(i),[-305 305],[-305 305]);    
    h=NS_AddFrameForArrayLayout(1,LineWidth);
    h=gca;
    set(h,'Visible','off');
end