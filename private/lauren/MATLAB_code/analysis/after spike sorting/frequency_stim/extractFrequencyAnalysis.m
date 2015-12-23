function [latencies successes stimAmp shortName] = extractFrequencyAnalysis(elecRespPath, movieInd, nSequence)

temp = load(elecRespPath);


stimAmp = temp.elecResp.stimInfo.stimAmps(movieInd);

%check that this analysis is finalized
if ~temp.elecResp.analysis.finalized(movieInd) && ~isempty(temp.elecResp.analysis.type{movieInd})
    warning(['analysis for ' elecRespPath ' movie ' num2str(temp.elecResp.stimInfo.movieNos(movieInd))  ' hasn''t been finalized'])
end



%latencies = zeros(5,5)
latencies = temp.elecResp.analysis.latencies{movieInd};

if exist('nSequence', 'var')
    latencies = reshape(latencies, nSequence, []); %converts to pulses in sequence x repetitions
end %if nSequence not specified, leave as total pulse number x 1
    
successes = latencies~=0;

end