clear;
NS_GlobalConstants=NS_GenerateGlobalConstants(61);
p='E:\2008-08-26-0\';
FileName='data005\data005000.bin';
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('E:\analysis\2008-08-26-0\data005\data005.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('E:\analysis\2008-08-26-0\data005\data005.neurons');


%p='D:\2008-08-27-4\';
%FileName='data004\data004000.bin';
%paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('D:\2008-08-27-4\analysis\data004\data004.params');
%neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\2008-08-27-4\analysis\data004\data004.neurons');

idList = neuronFile.getIDList();

NeuronID=558;
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings1=spikeTimes(1,1:min(1000,length(spikeTimes)))-7;
Channels=[1:64];
CenterChannel=44;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
Radius=1;
Channels=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';
%Channels=[44 2];
[RAWtraces,signal]=NS_AverageTraces([p FileName],Timings1-1,Channels,[-3 57],NS_GlobalConstants);
[x,y]=NS_EstimateSomaPosition(signal,Channels);

DataPath='E:\analysis\2008-08-26-0\data008_proba4';
%DataPath='D:\analysis\2008-08-27-4\data005';
PatternNumber=295;
MovieNumber=104;
Amplitudes=NS_AmplitudesForPattern(DataPath,Channels,PatternNumber,MovieNumber,NS_GlobalConstants);
[V]=NS_EstimatePotential(Amplitudes,Channels,x,y,20)
Grad=NS_EstimatePotentialGradient(Amplitudes,Channels,x,y,20)
[xs,ys]=NS_EstimateMassCenter(Amplitudes,Channels);