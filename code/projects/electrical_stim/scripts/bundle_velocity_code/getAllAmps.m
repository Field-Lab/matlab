function amps = getAllAmps(pathToAnalysisData, patternNo)
%%Returns list of stimulus amplitudes associated with particular data file and pattern number

movieNos = findMovieNos(pathToAnalysisData, patternNo); 
numMovies = length(movieNos);
amps = zeros(1, numMovies);
for i = 2:numMovies; %first one is extra so skip it
	[foo, ~, ~] = getStimAmps(pathToAnalysisData, 284, movieNos(i)); 
	amps(i) = foo; 
end
