function checkSystematicResids(pathToData, patternNos, neuronID, movieNosRange)

cd(pathToData)

nPatterns = length(patternNos);

temp = load(['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(1))]);
nMovies = length(temp.elecResp.stimInfo.movieNos);
nMoviesAll = length(movieNosRange);

residuals = zeros(nPatterns, nMoviesAll);
residualsUsedBin = false(nPatterns, nMoviesAll);
residualsMeans = zeros(1, nMoviesAll);
residualsStDevs = zeros(1, nMoviesAll);
amps = zeros(1, nMoviesAll);
ampDivisions = [];

for ii = 1:nPatterns
    ii
    temp = load(['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(ii))]);
    elecResp = temp.elecResp;
    a = elecResp.analysis.erfParams(1);
    b = elecResp.analysis.erfParams(2);
    
    %for jj = 1:nMovies
    for jj = 1:nMoviesAll
        movieInd = find(elecResp.stimInfo.movieNos == movieNosRange(jj));
        if ~isempty(movieInd)
            if ~isempty(elecResp.analysis.type{movieInd})
                x = abs(elecResp.stimInfo.stimAmps(movieInd));
                
                if amps(jj) ~= 0
                    if amps(jj) ~= x
                        disp(['warning--amplitudes for movie ' num2str(jj), ' not all equal'])
                    end
                else
                    amps(jj) = x;
                end 
                    
                %only include residuals that are in sensitive portion of curve (+/- 1 SD)
                if x >= elecResp.analysis.threshold*(1-1/(-sqrt(2)*b)) && x <= elecResp.analysis.threshold*(1+1/(-sqrt(2)*b))
                    expected = 0.5*(1 + erf(a*x+b));
                    residuals(ii,jj) = expected - elecResp.analysis.successRates(movieInd);
                    if residuals(ii,jj) == 0
                        keyboard
                    end
                    
                    residualsUsedBin(ii,jj) = true;
                end
            end
        end
    end
end

prevAmp = 0;
for jj = 1:nMoviesAll
    residualsUsed = residuals(residualsUsedBin(:,jj), jj);
    residualsMeans(jj) = mean(residualsUsed);
    residualsStDevs(jj) = std(residualsUsed);
    if amps(jj) ~= prevAmp && amps(jj) ~=0
        ampDivisions = [ampDivisions jj-0.5]; %#ok<AGROW>
        prevAmp = amps(jj);
    end
end

figure('position', [100 100 1000 250])
hold on
for jj = 1:nMoviesAll
    for ii = 1:nPatterns
        if residualsUsedBin(ii,jj)
            plot(jj, residuals(ii,jj), 'k.')
        end
    end
    plot(jj, residualsMeans(jj), 'ro')
    plot([jj jj], [residualsMeans(jj)-residualsStDevs(jj) residualsMeans(jj)+residualsStDevs(jj)], 'r-')
end

for ii = 1:length(ampDivisions)
    plot([ampDivisions(ii) ampDivisions(ii)], [-0.5 0.5], 'b-')
end

text(5, 0.7, [pathToData ' neuron ' num2str(neuronID)])
set(gca, 'ylim', [-0.8, 0.8])

xlabel(['movie number', 10, '(vertical lines indicate grouping by stimulus amplitude)'])
ylabel(['residuals (from erf fit),', 10, 'with means +/- SD in red'])

hold off
