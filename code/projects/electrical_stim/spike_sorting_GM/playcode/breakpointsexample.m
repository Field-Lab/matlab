 for p=1:58
patternNo=patterns(p);

movieNos           = findMovieNos([pathToAnalysisData],patternNo);

for m=1:length(movieNos)
    [bbb ch cc dd ee ff]=getStimAmps([pathToAnalysisData], patternNo, movieNos(m))
  
    ampi{p}(m,:)=ff;
end

 end