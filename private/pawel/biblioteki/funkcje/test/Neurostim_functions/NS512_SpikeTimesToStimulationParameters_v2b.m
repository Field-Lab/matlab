function SpikeTimesCombined=NS512_SpikeTimesToStimulationParameters_v2b(NeuronFilePath,MovieFilePath,DuplicatesFilePath,OutputPath);
%Ta funkcja generuje dla danego neuronu (powiedzmy o numerze 15) plik
%wyjsciowy o nazwie 'ID=15'. Plik ten zawiera macierz o wymiarach 5xN,
%gdzie N - ilo?c spik�w dla danego neuronu. UWAGA: m�wimy tutaj o fizycznym
%neuronie, na kt�ry moze sie skladac kilka neurons ID znalezione przez
%Vision. Stad tez w parametrach wejsciowych jest sciezka do pliku, kt�ry
%zawiera iinformacje o duplikatach (zidentyfikowanych przez Vision).
NS_GlobalConstants=NS_GenerateGlobalConstants(500);

%paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile(ParamsFilePath);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(NeuronFilePath);

PrimaryNeurons=NS512_IdentifyPrimaryNeurons_v2(NeuronFilePath,DuplicatesFilePath)
length(PrimaryNeurons)

fid=fopen(DuplicatesFilePath);
C=textscan(fid,'%d%s%s%d%s%f32');
DuplicateNeurons=C{1};
fclose(fid);

MoviesBegins=NS512_MoviesBegins(MovieFilePath,NS_GlobalConstants);

for neuron=1:length(PrimaryNeurons)
    neuron
    PrimaryNeuronID=PrimaryNeurons(neuron)
    DuplicatesOfGivenPrimaryNeuron=find(C{4}==PrimaryNeuronID)
    duplicateID=DuplicateNeurons(DuplicatesOfGivenPrimaryNeuron)'
    
    Duplicates=[PrimaryNeuronID duplicateID]
    %for i=1:sa(1) %unfortunately, the loop is required    
    %    PrimaryNeuron=str2num(a.textdata{i,4});
    %    if PrimaryNeuron==PrimaryNeuronID
    %        Duplicates=[Duplicates str2num(a.textdata{i,1})];
    %   end            
    %end
    
    spikeTimes=[];
    for j=1:length(Duplicates)
        spikeTimes = [spikeTimes double(neuronFile.getSpikeTimes(PrimaryNeuronID))']; % for given neuron, import the spikes times
    end
    SortedSpikeTimes=sort(unique(spikeTimes));
    SpikeToSpikeTime=diff(SortedSpikeTimes);    
    for s=1:length(SpikeToSpikeTime)        
        if SpikeToSpikeTime(s)<10
            SortedSpikeTimes(s)=0;
        end
    end
    length(SortedSpikeTimes)
    length(find(SortedSpikeTimes>0))
    SpikeTimesCombined=unique(SortedSpikeTimes(find(SortedSpikeTimes>0)));
    find(SpikeTimesCombined==0)
    L=length(spikeTimes);
                
    dane=zeros(5,length(SpikeTimesCombined));
    for i=1:length(SpikeTimesCombined)        
        [dane(1,i),dane(2,i),dane(3,i),dane(4,i),dane(5,i),dane(6,i)]=NS512_SpikeTimesToStimulationPatterns_v3b(MovieFilePath,SpikeTimesCombined(i),MoviesBegins,NS_GlobalConstants);
    end
    FullName=[OutputPath '\ID=' num2str(PrimaryNeuronID) 'c'];
    fid=fopen(FullName,'wb','ieee-le');                                    
    fwrite(fid,dane,'int32');
    fclose(fid);
end

% to read:
%fid=fopen('ID=3','r','ieee-le'); 
%a=fread(f,'int32');
%l=length(a);
%b=reshape(a,4,l/4);