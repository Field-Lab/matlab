function [spikeTimesToUse order] = spike_train_median_calculator(spikeTimesAll, nToChoose, choiceStart, costParam, normType, trainLength, varargin)


p = inputParser;

p.addRequired('spikeTimesAll', @iscell)
p.addRequired('nToChoose', @isnumeric)
p.addRequired('choiceStart', @isnumeric)
p.addRequired('costParam', @isnumeric)
p.addRequired('normType', @ischar)
p.addRequired('trainLength', @isnumeric) %for plotting

p.addParamValue('cellIDs', [], @isnumeric) %used for labeling plots only


p.parse(spikeTimesAll, nToChoose, choiceStart, costParam, normType, trainLength, varargin{:})

cellIDs = p.Results.cellIDs;



%% basics

nCells = length(spikeTimesAll{1});
nReps = length(spikeTimesAll);

nSpikes = zeros(nReps, nCells);
minISI = 1;
for ii = 1:nReps
    for jj = 1:nCells
        nSpikes(ii,jj) = length(spikeTimesAll{ii}{jj});
        minISI = min([minISI min(diff(spikeTimesAll{ii}{jj}))]);
    end
end
nSpikeMeans = mean(nSpikes, 1); %mean number of spikes for each cell


%% sniff tests (use if messing with algorithm or parameters)

if 0    %add spikes at same time to all repetitions of originals
    for i = 1:nCells
        randSpikes = 0.8*rand(5+i,1);
        for k = 1:nReps
            spikeTimesAll{k}{i} = [spikeTimesAll{k}{i}; randSpikes];
        end
    end
end

if 0 % replace one cell's spike train with single spikes at random times and leave the rest of the cells empty
    i = 1;
    spikeTimesCell = [];
    for k = 1:nReps
        spikeTimesAll{k}{i} = 0.08*rand(1);
        spikeTimesCell = [spikeTimesCell spikeTimesAll{k}{i}];
    end
    for i = 2:nCells
        for k = 1:nReps
            spikeTimesAll{k}{i} = [];
        end
    end
end

if 0 %replace one cell's spike train with bursts of different numbers of spikes at roughly the same time
    i = 1;
    spikeCountCell = [];
    for k = 1:nReps
        randVal = 0.24*rand(1);

        %spikeTimesAll{k}{i} = 0.5:0.01:(0.5 + randVal*2);
        
        spikeTimesAll{k}{i} = (0.5-randVal+0.02*rand(1)):0.02:(0.5+randVal);
        spikeTimesAll{k}{i} = spikeTimesAll{k}{i} + 0.004*rand(1);
        spikeCountCell = [spikeCountCell length(spikeTimesAll{k}{i})];
    end
    for i = 2:nCells
        for k = 1:nReps
            spikeTimesAll{k}{i} = [];
        end
    end
end

%% check against max likelihood estimate (alternative measure of spike
% train similarity)

if 0
    kernWidth = 0.01; %SD of Guassian kernel used to estimate continuous spike rate, in seconds
    
    [logLikes spikeRates times] = max_likelihood_train_calculator(spikeTimesAll, kernWidth);
    
    % figure
    % for ii = 1:nCells
    %     hold on
    %     subplot(nCells,1,ii)
    %     plot(times, spikeRates{ii})
    %     hold off
    % end
    
    [tmp MLorder] = sort(logLikes, 2, 'descend');
        
    % % slider plot
    % sliderFig = figure;
    % slider = make_loop_slider_list(1,1,nReps);
    %
    % while ishandle(sliderFig)
    %     ii = round(get(slider,'Value'));
    %     for jj = 1:nCells
    %         subplot(nCells+1,1,jj+1)
    %         cla
    %         hold on
    %         plot(times, spikeRates{jj})
    %         for kk = 1:length(spikeTimesAll{ii}{jj})
    %             plot([spikeTimesAll{ii}{jj}(kk) spikeTimesAll{ii}{jj}(kk)], [0 0.5], 'k-')
    %         end
    %         xlabel('time (s)')
    %         ylabel('estimated spike rate with raster')
    %         hold off
    %     end
    %
    %     subplot(nCells+1,1,1)
    %     cla
    %     hold on
    %     plot(logLikes, 'k-')
    %     plot(ii,logLikes(ii), 'ro')
    %     xlabel('trial')
    %     ylabel('log likelihood')
    %     hold off
    %
    %     uiwait;
    % end
end

%% compute pair-wise Victor spike train distances (costs)


costParam = costParam/1000; %in seconds; the distance that a spike moves that is equivalent (in cost) to deleting and inserting spike

costParam = 2/costParam; %cost per second to move spike


costs = zeros(nReps, nReps);
sumCost = zeros(1,nReps);
for i = 1:nReps %trial being compared
    for j = 1:nReps %trial being compared TO
        pairCost = 0;
        for k = 1:nCells %calculate pairwise cost for each cell and sum
            d = spkd(spikeTimesAll{i}{k}, spikeTimesAll{j}{k}, costParam);

            if strcmpi(normType, 'none')
                thisCost = d; %no normalization
            elseif strcmpi(normType, 'origin')
                thisCost = d/length(spikeTimesAll{i}{k}); %normalize by number of spikes in spike train being compared
            elseif strcmpi(normType, 'destination')
                thisCost = d/length(spikeTimesAll{j}{k}); %normalize by number of spikes in spike train being compared TO
            elseif strcmpi(normType, 'both')
                thisCost = d/(length(spikeTimesAll{i}{k}) + length(spikeTimesAll{j}{k})); %normalize by total number of spikes in both spike trains
            else
                error('unrecognized normalization type')
            end
            
            sumCost(i) = sumCost(i) + thisCost; %sum over distances of all trains from ith train
            pairCost = pairCost + thisCost; %sum over cells for ith vs. jth trains
            
        end
        costs(i,j) = pairCost;
    end
end

meanCosts = sumCost;


%% plotting options

if 0 %plots histogram of mean distances
    figure
    hist(meanCosts)
end


if 0 %plots linkage map
    %converts matrices to "pdist" format
    %note that this is only valid if costs matrix is symmetric (if cost between a pair of spike
    %trains is normalized that same way in both "directions")
    dists = [];
    for i = 1:nReps
        dists = [dists costs(i,i+1:end)];
    end

    linkmap = linkage(dists);

    figure
    dendrogram(linkmap, 0);
end


%% choose 'median' spike trains

[y order] = sort(meanCosts);

spikeTimesToUse = cell(nToChoose, 1);
for ii = 1:nToChoose
    spikeTimesToUse{ii} = spikeTimesAll{order(ii+choiceStart-1)};
end

%% compare ordering based on distance metric vs. likelihood

% figure
% plot(meanCosts, logLikes, 'ko')
% xlabel(['mean distance, normalization type = ' normType])
% ylabel('log likelihood')

%% compare ordering to number of spikes in trial

figure
plot(meanCosts, sum(nSpikes,2), 'ko')
xlabel(['mean distance, normalization type = ' normType])
ylabel('number of spikes in trial')

% figure
% plot(logLikes, sum(nSpikes,2), 'ko')
% xlabel('log likelihood')
% ylabel('number of spikes in trial')


%% plot rasters, ordered according to each metric

if 1
    figure
    for i = 1:nCells
        axes('position', [0.1 (nCells - i + 1)/(nCells+1) 0.8 1/(nCells+2)])
        hold on
        for k = 1:nReps
            for j = 1:length(spikeTimesAll{order(k)}{i})
                plot([spikeTimesAll{order(k)}{i}(j) spikeTimesAll{order(k)}{i}(j)], [k-1 k], 'k-', 'LineWidth', 1)
            end
        end
        hold off

        set(gca, 'yLim', [0 nReps], 'xlim', [0 trainLength])
        if i == nCells
            xlabel('time (s)')
        else
            set(gca, 'xtick', [])
        end
        if ~isempty(cellIDs)
            ylabel(['cell' num2str(cellIDs(i))])
        end
        title(['trial order, based on distance metric (bottom to top), normalization type = ' normType])
    end
    
%     figure
%     for i = 1:nCells
%         axes('position', [0.1 (nCells - i + 1)/(nCells+1) 0.8 1/(nCells+2)])
%         hold on
%         for k = 1:nReps
%             for j = 1:length(spikeTimesAll{MLorder(k)}{i})
%                 plot([spikeTimesAll{MLorder(k)}{i}(j) spikeTimesAll{MLorder(k)}{i}(j)], [k-1 k], 'k-', 'LineWidth', 1)
%             end
%         end
%         hold off
%         
%         set(gca, 'yLim', [0 nReps], 'xlim', [0 trainLength])
%         if i == nCells
%             xlabel('time (s)')
%         else
%             set(gca, 'xtick', [])
%         end
%         if ~isempty(cellIDs)
%             ylabel(['cell' num2str(cellIDs(i))])
%         end
%         title('trial order, based on log likelihood (bottom to top)')
%     end
end


