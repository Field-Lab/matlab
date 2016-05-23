NeuronFilePath = 'D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002\data002.neurons';
DuplicatesFilePath='D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_duplicates.txt';
PrimaryNeurons=NS512_IdentifyPrimaryNeurons_v2(NeuronFilePath,DuplicatesFilePath);


SpikeTimesPath='C:\home\Pawel\nauka\512stim_paper\NeuronAnalysis\2010-09-14-0\data002\SpikeTimesCombined';
GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\NeuronAnalysis\2010-09-14-0\data002\GaussesFiles';

SPFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data002\sp_files_2015-05-24\';
GaussesFilesPath='C:\pawel\nauka\512paper\SpikesAnalysis\GaussesFiles2\';



TimeOffset=20;
MaxFitError=0.1;
%[CorrData]=NS512_SpikeCorrelationsMultielectrodeStimulation(PrimaryNeurons,SpikeTimesPath,GaussesFilesPath,TimeOffset,MaxFitError);

s1=sign(CorrData(:,5)-1*CorrData(:,6));
s2=sign(CorrData(:,5)-1*CorrData(:,7));

s3=s1+s2;
length(find(s3==2))