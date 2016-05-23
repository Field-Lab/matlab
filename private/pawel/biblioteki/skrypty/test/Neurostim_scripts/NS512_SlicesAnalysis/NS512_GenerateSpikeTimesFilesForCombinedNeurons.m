%SpikeTimesCombined=NS512_SpikeTimesToStimulationParameters_v2('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001\data001.neurons','I:\analysis\slices\2013-12-12-3-PH\movie001','D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\duplicates2.txt','D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\SpikeTimesCombined');

NeuronFilePath = 'D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002\data002.neurons';
MovieFilePath='I:\backup1\dane\2010-09-14-0\movie002';
DuplicatesFilePath='D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_duplicates.txt';
OutputPath='D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_Matlab\SpikeTimesCombined';

SpikeTimesCombined=NS512_SpikeTimesToStimulationParameters_v2(NeuronFilePath,MovieFilePath,DuplicatesFilePath,OutputPath);

%NS512_RasterPlots_1electrode_scan_2014_05_20;
