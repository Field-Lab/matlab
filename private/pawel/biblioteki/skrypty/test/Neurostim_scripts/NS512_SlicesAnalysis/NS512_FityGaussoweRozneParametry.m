% Spikes analysis, szczur
zwierze=1;
NumberOfAmplitudes=17;

SPFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data002\sp_files-2015-05-24\';
MovieFilePath='C:\home\Pawel\nauka\512stim_paper\MovieFiles\2010-09-14-0\movie002';

FiguresPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data002\GaussesFigures-2015-06-11\';
GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data002\GaussesFiles-2015-06-11\';
WersjaFitowania=1;
test1;

FiguresPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data002\GaussesFigures-2015-06-11b\';
GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data002\GaussesFiles-2015-06-11b\';
WersjaFitowania=2;
test1;

SPFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data005\sp_files-2015-05-24\';
MovieFilePath='C:\home\Pawel\nauka\512stim_paper\MovieFiles\2010-09-14-0\movie005';

FiguresPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data005\GaussesFigures-2015-06-11\';
GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data005\GaussesFiles-2015-06-11\';
WersjaFitowania=1;
test1;

FiguresPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data005\GaussesFigures-2015-06-11b\';
GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data005\GaussesFiles-2015-06-11b\';
WersjaFitowania=2;
test1;

% Spikes analysis, mysz
clear
zwierze=2;
NumberOfAmplitudes=25; 

MovieFilePath='C:\home\Pawel\nauka\512stim_paper\MovieFiles\2013-12-12-3\movie001';
SPFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\sp_files-2015-05-24\';

FiguresPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\GaussesFigures-2015-06-11\';
GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\GaussesFiles-2015-06-11\'; % jedyna roznica w porownaniu do 2015-05-24: 
WersjaFitowania=1;
test1;

FiguresPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\GaussesFigures-2015-06-11b\';
GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\GaussesFiles-2015-06-11b\'; % jedyna roznica w porownaniu do 2015-05-24: 
WersjaFitowania=2;
test1;

% Neuron analysis, szczur, data002 only
break
clear
GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\NeuronAnalysis\2010-09-14-0\data002\GaussesFiles-2015-06-11\';
FiguresPath='C:\home\Pawel\nauka\512stim_paper\NeuronAnalysis\2010-09-14-0\data002\GaussesFigures-2015-06-11\';
WersjaFitowania=1;
Stim512Paper_NeuronAnalysis.m

GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\NeuronAnalysis\2010-09-14-0\data002\GaussesFiles-2015-06-11b\';
FiguresPath='C:\home\Pawel\nauka\512stim_paper\NeuronAnalysis\2010-09-14-0\data002\GaussesFigures-2015-06-11b\';
WersjaFitowania=2;
Stim512Paper_NeuronAnalysis.m