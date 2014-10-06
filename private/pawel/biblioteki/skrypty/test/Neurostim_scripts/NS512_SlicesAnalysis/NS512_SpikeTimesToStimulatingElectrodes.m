NS_GlobalConstants=NS_GenerateGlobalConstants(500);

paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\Vision_output\2010-09-14-0\data002\data002.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\Vision_output\2010-09-14-0\data002\data002.neurons');
idList = neuronFile.getIDList(); %this imports list of neurons IDs from the file that is defined above as 'neuronFile';

full_path='D:\Home\Data\slices\2010-09-14-0\movie002';
MoviesBegins=NS512_MoviesBegins(full_path,NS_GlobalConstants);

for neuron=1:length(idList)
    neuron
    NeuronID=idList(neuron);
    spikeTimes = double(neuronFile.getSpikeTimes(NeuronID))'; % for given neuron, import the spikes times
    dane=zeros(4,length(spikeTimes));
    for i=1:length(spikeTimes)
        %i
        [dane(1,i),dane(2,i),dane(3,i),dane(4,i)]=NS512_SpikeTimesToStimulationPatterns(full_path,spikeTimes(i),MoviesBegins,NS_GlobalConstants);
    end
    FullName=['C:\pawel\nauka\analiza\slices\2010-09-14-0\data002minus009paramSninya\2010-09-14-0\dane\ID=' num2str(NeuronID)];
    FullName=['C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\dane\ID=' num2str(NeuronID)];
    fid=fopen(FullName,'wb','ieee-le');                                    
    fwrite(fid,dane,'int32');
    fclose(fid);
end

% to read:
%fid=fopen('ID=3','r','ieee-le'); 
%a=fread(f,'int32');
%l=length(a);
%b=reshape(a,4,l/4);