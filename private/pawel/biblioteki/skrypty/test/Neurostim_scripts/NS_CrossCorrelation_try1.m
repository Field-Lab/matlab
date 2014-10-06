FileName='005000.bin';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\analysis\2008-08-26-0\data005\data005.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\analysis\2008-08-26-0\data005\data005.neurons');
idList = neuronFile.getIDList();
NeuronID=376;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings1=spikeTimes(1,1:min(1000,length(spikeTimes)))-35;
Channels=26;
cd D:\2008-08-26-0\data005;
[RAWtraces,signal]=NS_AverageTraces(FileName,Timings1-1,Channels,[-2 157],NS_GlobalConstants);
figure(1)
plot(signal)

SC=zeros(1000,319);

for i=1:1000
    s1=RAWtraces(i,:);
    s1=s1-mean(s1);
    s2=signal-mean(signal);
    SC(i,:)=xcorr(s1,s2);
end

figure(1)
plot(RAWtraces(1:200,:)');
figure(2)
plot(SC(1:200,:)')