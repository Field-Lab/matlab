function [likelihoods spikeRates times] = max_likelihood_train_calculator(spikeTimes, kernWidth)

binSize = 0.01; %in s


nTrials = length(spikeTimes);
nCells = length(spikeTimes{1});

maxSpikeTime = 0;

for ii = 1:nTrials
    for jj = 1:nCells
        maxSpikeTime = max([maxSpikeTime max(spikeTimes{ii}{jj})]);
    end
end

%continuous spike rate estimate = spike trains convolved with Gaussian kernels = sum
%of normalized Gaussians centered at spike times divided by number of
%trials


%sampleSize = 0.0001; %in seconds
%times = 0:sampleSize:maxSpikeTime+sampleSize;
times = binSize:binSize:maxSpikeTime+binSize; %upper edges of bins

nBins = ceil(max(times)/binSize);

spikeRates = cell(nCells,1);

for jj = 1:nCells
    spikeRates{jj} = zeros(1, nBins);
end


% for jj = 1:nCells
%     spikeRates{jj} = zeros(1, length(times));
% end

%figure
for ii = 1:nTrials
    for jj = 1:nCells
        %subplot(2,1,jj)
        %hold on
        
        for kk = 1:nBins
            binStart = (kk-1)*binSize;
            binEnd = kk*binSize;
            if any(spikeTimes{ii}{jj} > binStart & spikeTimes{ii}{jj} <= binEnd)
                spikeRates{jj}(kk) = spikeRates{jj}(kk) + 1/nTrials;
            end

        end
        
%         for kk = 1:length(spikeTimes{ii}{jj})%go through each spike (probably much more efficient way to do this)
%             spikeRates{jj} = spikeRates{jj} + (normpdf(times, spikeTimes{ii}{jj}(kk), kernWidth))/nTrials;
%             %plot(spikeRates{jj})
%         end
        %hold off
    end
end

%hold off

%% determine likelihood of each observed trial

likelihoods = zeros(1,nTrials);

for ii = 1:nTrials
    for jj = 1:nCells
        
        %         for kk = 1:length(spikeTimes{ii}{jj})
        %
        %             %determine which time bin spike falls in
        %             iBin = round(spikeTimes{ii}{jj}(kk)/sampleSize)+1;
        %
        %             likelihoods(ii) = likelihoods(ii) + log(spikeRates{jj}(iBin));
        %
        %         end
        
        for kk = 1:nBins
            binStart = (kk-1)*binSize;
            binEnd = kk*binSize;
            if any(spikeTimes{ii}{jj} > binStart & spikeTimes{ii}{jj} <= binEnd)
                likelihoods(ii) = likelihoods(ii) + log(spikeRates{jj}(kk));
            else
                likelihoods(ii) = likelihoods(ii) + log(1-spikeRates{jj}(kk));
            end
        end
    end
end


