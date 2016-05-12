%Parameters
fname = '/Volumes/Analysis/2015-10-06-3/data001-data002/'
patternNo = 284;

%Load EI after stimulation
[rawData, amps] = generateEiFromStimPattern(fname, patternNo, 'suppressPlots', true);

%Get Stim Amps
amps = getAllAmps(fname, patternNo);
namps = length(amps) - 1; %function returns one more than exists for some reason

%Get involved electrodes at each amp
eleclist = cell(1, namps);
for i = 1:namps
	elecs = getBundleElecs(rawData, i, -50);
	eleclist{i} = elecs;
end

