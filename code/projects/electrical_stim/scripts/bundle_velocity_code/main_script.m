%Parameters
fname = '/Volumes/Analysis/2015-10-06-3/data001-data002/'
patternNo = 284;

%Options
plflg = 0; 

%Load EI after stimulation
if ~exist('rawData', 'var')
	[rawData, amps] = generateEiFromStimPattern(fname, patternNo, 'suppressPlots', true);
	%Subtract global mean
	rawData = rawData - mean2(rawData);
end


%Get Stim Amps
amps = getAllAmps(fname, patternNo);
namps = length(amps) - 1; %function returns one more than exists for some reason

%Get involved electrodes at each amp
if plflg == 1;
	eleclist = cell(1, namps);
	for i = 1:namps
		elecs = getBundleElecs(rawData, i, -50);
		eleclist{i} = elecs;
	end

	%plot electrodes
	plotElecsDots(patternNo, eleclist{30})
end

%Get involved electrodes at specific time and find clusters
elecs = getBundleElecs_overtime(rawData, 38, -50, 23);
plotElecsDots(patternNo, elecs)
howManyClusters(elecs, patternNo)
