NS_GlobalConstants=NS_GenerateGlobalConstants(500);

paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';
idList=[96 168 197 216 217 261];
idList=[512 513 529 558 573]
idList=[96 216 217 289 183 261 800 858 1022 1097];

idList=289

full_path='I:\analysis\slices\2013-12-12-3-PH\movie005';
MoviesBegins=NS512_MoviesBegins(full_path,NS_GlobalConstants);

for neuron=1:length(idList)
    neuron
    NeuronID=idList(neuron);
    spikeTimes = double(neuronFile.getSpikeTimes(NeuronID))'; % for given neuron, import the spikes times
    dane=zeros(5,length(spikeTimes));
    for i=1:length(spikeTimes)
        %i
        [dane(1,i),dane(2,i),dane(3,i),dane(4,i),dane(5,i)]=NS512_SpikeTimesToStimulationPatterns_v3(full_path,spikeTimes(i),MoviesBegins,NS_GlobalConstants);
    end
    %FullName=['C:\pawel\nauka\analiza\slices\2010-09-14-0\data002minus009paramSninya\2010-09-14-0\dane\ID=' num2str(NeuronID)];
    FullName=['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\Matlab\dane\ID=' num2str(NeuronID)];
    fid=fopen(FullName,'wb','ieee-le');                                    
    fwrite(fid,dane,'int32');
    fclose(fid);
end

% to read:
%fid=fopen('ID=3','r','ieee-le'); 
%a=fread(f,'int32');
%l=length(a);
%b=reshape(a,4,l/4);