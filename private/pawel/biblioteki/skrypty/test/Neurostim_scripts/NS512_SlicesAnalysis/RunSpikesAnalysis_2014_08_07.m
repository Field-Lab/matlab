clear
MovieFile='G:\data\2013-12-12-3-PH\movie001';
DataPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc';
ArtifactDataPath=DataPath;
FigurePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\SpikeAnalysis\2012-12-12-3-PH-2014-08-07';

NumberOfAmplitudes=25; %25 for 2013-12-12-3
NumberOfMovieSequences=16; %16 for 2013-12-12-3
NumberOfMovieSequences2=18; %18 for 2013-12-12-3

SzybkiTestDetekcjiSpikow_v7.m;
NumberOfStimulatedElectrodesMouse=NumberOfStimulatedElectrodes;
StimulatedElectrodesMouse=StimulatedElectrodes;
save C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\StimulatedElectrodesMouse StimulatedElectrodesMouse;

break
clear

MovieFile='I:\analysis\slices\2010-09-14-0\TTX_sub\movie002';
DataPath='G:\analysis\2010-09-14-0\data002_preproc';
ArtifactDataPath=DataPath;
FigurePath='D:\Home\Pawel\analysis\slices\2010-09-14-0\SpikesAnalysis\data002';

NumberOfAmplitudes=17; %25 for 2013-12-12-3
NumberOfMovieSequences=8; %16 for 2013-12-12-3
NumberOfMovieSequences2=8; %18 for 2013-12-12-3

SzybkiTestDetekcjiSpikow_v7.m;
NumberOfStimulatedElectrodesRat=NumberOfStimulatedElectrodes;
save C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\StimulatedElectrodesRat StimulatedElectrodesRat;