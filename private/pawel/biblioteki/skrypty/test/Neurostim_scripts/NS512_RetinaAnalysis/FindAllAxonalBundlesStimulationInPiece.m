Patterns=[1:512];
Eletrodes=[1:512];
Movies=[1:2:63];
Threshold=20;
SamplesToAnalyze=[8:60];

tic
DataPath='G:\D\Home\Pawel\analysis\retina\2012sept\2012-09-27-4\scan_proc';
Thresholds=NS512_FindAllAxonalBundleStimulation(DataPath,Patterns,Movies,Threshold,SamplesToAnalyze);
cd G:\D\Home\Pawel\analysis\retina\2012sept\2012-09-27-4;
save Thresholds

DataPath='G:\D\Home\Pawel\analysis\retina\2012sept\2012-09-24-0\scan_proc';
Thresholds=NS512_FindAllAxonalBundleStimulation(DataPath,Patterns,Movies,Threshold,SamplesToAnalyze);
cd G:\D\Home\Pawel\analysis\retina\2012sept\2012-09-24-0;
save Thresholds

DataPath='G:\D\Home\Pawel\analysis\retina\2012sept\2012-09-24-3\scan_proc';
Thresholds=NS512_FindAllAxonalBundleStimulation(DataPath,Patterns,Movies,Threshold,SamplesToAnalyze);
cd G:\D\Home\Pawel\analysis\retina\2012sept\2012-09-24-3;
save Thresholds
toc