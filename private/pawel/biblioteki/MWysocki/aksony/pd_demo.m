InputPath = 'D:\Home\Pawel\analysis\2012_09\2012-09-27-4\scan_new';
OutputPath = 'D:\Home\Pawel\analysis\retina\2012-09-27-4\analysis_2014_05_16\MWysocki';
DataID = 'retina1';
PatternRange = [1:512];
MovieRange = [63:-2:51];
MovieNumber = 63;
BadElectrodes = [420 353 354 488];
%BadElectrodes = getBadElectrodes( InputPath, PatternRange, MovieNumber)
InitialDirection = [];

file = PD_AcquireData(InputPath, OutputPath, DataID, PatternRange, MovieRange, BadElectrodes, InitialDirection);