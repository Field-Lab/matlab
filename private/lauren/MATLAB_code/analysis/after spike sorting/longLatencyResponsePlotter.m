clear all

[neuron special_neuron] = cell_list_low_freq_stim();

nNormalNeurons = length(neuron);

% raster plots

histRateMin = 0.5; %response rate (from curve fit) above which long-latency responses are collected for histogram plot
histRateMax = 0.9;

for ii = 1:length(neuron)
    datarun.names.rrs_neurons_path = neuron(ii).path;
    datarun = load_neurons(datarun);
    
    load(neuron(ii).movieTimePath);
    
    nMovies = sum(diff(movieNumberTimes(1,:))~=0)+1;
    
    load(neuron(ii).elecRespPath)
    elecResp = checkForUnfinishedAnalysis(elecResp, 100);
    
    thresh = elecResp.analysis.threshold;
    sig = 1/(elecResp.analysis.erfParams(1)*sqrt(2));

    %amplitude corresponding to histRateMin in cumulative Gaussian fit
    histAmpMin = norminv(histRateMin, thresh, sig);
    histAmpMax = norminv(histRateMax, thresh, sig);
    includedInHist = false(1,size(movieNumberTimes, 2));
        
    spikeTimesTemp = datarun.spikes{get_cell_indices(datarun, neuron(ii).id)}*20000; %in samples
    spikeTimes = cell(1,size(movieNumberTimes, 2));
    spikeTimesCat = []; %collect spike times for movies corresponding to above threshold amplitudes for histogram
    for kk = 1:size(movieNumberTimes, 2)
        t = movieNumberTimes(2,kk);
        if kk < size(movieNumberTimes, 2)
            tNext = movieNumberTimes(2,kk+1);
            spikeTimes{kk} = spikeTimesTemp(spikeTimesTemp >= t & spikeTimesTemp < tNext);
        else
            spikeTimes{kk} = spikeTimesTemp(spikeTimesTemp >= t);
        end
        spikeTimes{kk} = (spikeTimes{kk} - t)/20; %zero to stimulus and convert to ms
        
        movieInd = find(elecResp.stimInfo.movieNos == movieNumberTimes(1,kk));
        if abs(elecResp.stimInfo.stimAmps(movieInd)) >= histAmpMin &&...
                abs(elecResp.stimInfo.stimAmps(movieInd)) <= histAmpMax
            spikeTimesCat = [spikeTimesCat; spikeTimes{kk}];
            includedInHist(kk) = true;
        end
    end
    
    % extract elecResp analysis
    if length(elecResp.stimInfo.movieNos) ~= nMovies
        error('different number of movies in elecResp than value listed in movieNumberTimes file')
    end
    
    succRates = elecResp.analysis.successRates;
    
    succColors = jet(101);
    
    figure; hold on
    for kk = 1:length(spikeTimes)
        for mm = 1:length(spikeTimes{kk})
%             if includedInHist(kk)
%                 plot([spikeTimes{kk}(mm) spikeTimes{kk}(mm)], 0.1*[kk-1 kk], 'k-', 'LineWidth', 2)
%             else
%                 plot([spikeTimes{kk}(mm) spikeTimes{kk}(mm)], 0.1*[kk-1 kk], 'k-', 'LineWidth', 1)
%             end
            if includedInHist(kk)
                plot(spikeTimes{kk}(mm), 0.1*kk, 'k.', 'MarkerSize', 10)
            else
                plot(spikeTimes{kk}(mm), 0.1*kk, 'k.', 'MarkerSize', 2)
            end
        end
        %plot color patch representive response rate
        movieInd = find(elecResp.stimInfo.movieNos == movieNumberTimes(1,kk));
        patch([-5 0 0 -5], 0.1*[kk-1 kk-1 kk kk], succColors(round(100*succRates(movieInd))+1,:), 'edgeColor', 'none')
        
        if kk < length(spikeTimes)
            nextMovieInd = find(elecResp.stimInfo.movieNos == movieNumberTimes(1,kk+1));
            if abs(elecResp.stimInfo.stimAmps(movieInd)) < 2*thresh && abs(elecResp.stimInfo.stimAmps(nextMovieInd)) >= 2*thresh
                plot([0 100], 0.1*[kk kk], 'k-')
            end
            
            if abs(elecResp.stimInfo.stimAmps(movieInd)) < 1.5*thresh && abs(elecResp.stimInfo.stimAmps(nextMovieInd)) >= 1.5*thresh
                plot([0 100], 0.1*[kk kk], 'b-')
            end
        end
        clear movieInd nextMovieInd
    end
    
    disp(['stimulated through ' num2str(max(abs(elecResp.stimInfo.stimAmps))/thresh) 'x threshold'])
    disp(['(response probability = ' num2str(normcdf(max(abs(elecResp.stimInfo.stimAmps)), thresh, sig))  ')'])
    
    binEdges = 0:5:100;
    %latenciesMs = spikeTimesCat/20;
    histVals = psthPlotterBase(gca, binEdges, spikeTimesCat);
    title(['n' num2str(neuron(ii).id) '(' neuron(ii).type ')'])
    set(gca, 'xlim', [-5 100])
    
    %plot acf from WN datarun
    
    datarunWN.names.rrs_neurons_path = neuron(ii).WN.neuronsPath;
    datarunWN = load_neurons(datarunWN);
    
    [times acf] = get_correlation(datarunWN, neuron(ii).WN.id);
    acf = fliplr(acf);
    acf = [0 acf];
    times = [0 times];
    
    plot(times*1000, acf*2, 'r-')
    
    pause
    
    %%% neuron with doublet responses
%     if neuron(ii).id == 694
%         nSpikes = zeros(1, length(spikeTimes));
%         for jj = 1:length(spikeTimes)
%             nSpikes(jj) = sum(spikeTimes{jj}>50);
%         end
%         keyboard
%     end
    
    close(gcf)
    
    clear elecResp movieNumberTimes datarun datarunWN
end


%%
%off midget examples that required gui analysis (elecResp) to analyze
%long-latency responses due to low spontaneous firing rates

%note that this data organization was hacked together very quickly and
%without much planning, and is only intended for checking a couple cells
%quickly!  It should be correct in its current state but relies on several
%assumptions about how the data is split into elecResp movies and trials


for ii = 1:length(special_neuron)
    
    %get info from original elecResp file
    load(special_neuron(ii).originalElecRespPath)
    oElecResp = elecResp; clear elecResp
    oElecResp = checkForUnfinishedAnalysis(oElecResp, 100);
    
    thresh = oElecResp.analysis.threshold;
    sig = 1/(oElecResp.analysis.erfParams(1)*sqrt(2));
    
    %amplitude corresponding to histRateMin in cumulative Gaussian fit
    histAmpMin = norminv(histRateMin, thresh, sig);
    histAmpMax = norminv(histRateMax, thresh, sig);

    
    files = dir(special_neuron(ii).splitElecRespPath);
    for jj = length(files):-1:1
        if isempty(strfind(files(jj).name, ['elecResp_n' num2str(special_neuron(ii).id) '_p' num2str(special_neuron(ii).stimElec)]))
            files(jj) = [];
        end
    end
    
    load(special_neuron(ii).repStartTimesPath)
    
    load(special_neuron(ii).movieTimePath)
    spikeTimes = cell(1,size(movieNumberTimes, 2));

    oMovieRep = 0; %keeps track of which repitition of original movie we're on for this amplitude
    %collect spike times from elecResp files
    for jj = 1:length(files)
        load([special_neuron(ii).splitElecRespPath filesep files(jj).name])
        elecResp = checkForUnfinishedAnalysis(elecResp, 100);
        movieNos = elecResp.stimInfo.movieNos;
        
        %om = num2str(files(jj).name(strfind(files(jj).name, 'om')+2:strfind(files(jj).name, '.mat')-1));

        for kk = 1:length(movieNos)
            
            latencies = elecResp.analysis.latencies{kk};
            for mm = 1:length(latencies) %loop through trials
                if repOutStartTimes(movieNos(kk),mm) == 0
                    oMovieRep = oMovieRep+1; %advance movie rep number when value in repOutStartTimes goes back to 0, signifying first time interval for new movie rep
                    spikeTimes{oMovieRep} = [];
                end
                if latencies(mm)~=0
                    spikeTimes{oMovieRep} = [spikeTimes{oMovieRep} (latencies(mm)+repOutStartTimes(movieNos(kk),mm))/20];
                end 
            end
        end
        clear elecResp
    end
    

    %go back through and collect spike times from stim amplitudes that
    %satisfy histRateMin criterion
    includedInHist = false(1,length(spikeTimes));

    spikeTimesCat = [];
    for jj = 1:length(spikeTimes)
        spikeTimes{jj}(diff(spikeTimes{jj})<0.2) = []; %remove spikes that show up twice because they exist in overlap regions
        movieInd = find(oElecResp.stimInfo.movieNos == movieNumberTimes(1,jj));
        if abs(oElecResp.stimInfo.stimAmps(movieInd)) >= histAmpMin &&...
                abs(elecResp.stimInfo.stimAmps(movieInd)) <= histAmpMax
            spikeTimesCat = [spikeTimesCat; spikeTimes{jj}];
            includedInHist(jj) = true;
        end
    end
    
    % plot
    
    succRates = oElecResp.analysis.successRates;
    succColors = jet(101);
    
    figure; hold on
    for kk = 1:length(spikeTimes)
        for mm = 1:length(spikeTimes{kk})
            %plot([spikeTimes{kk}(mm) spikeTimes{kk}(mm)], 0.1*[kk-1 kk], 'k-', 'LineWidth', 1)
            if includedInHist(kk)
                plot(spikeTimes{kk}(mm), 0.1*kk, 'k.', 'MarkerSize', 10)
            else
                plot(spikeTimes{kk}(mm), 0.1*kk, 'k.', 'MarkerSize', 2)
            end
        end

        %plot color patch representive response rate
        movieInd = find(oElecResp.stimInfo.movieNos == movieNumberTimes(1,kk));
        patch([-5 0 0 -5], 0.1*[kk-1 kk-1 kk kk], succColors(round(100*succRates(movieInd))+1,:), 'edgeColor', 'none')
        
        if kk < length(spikeTimes)
            nextMovieInd = find(oElecResp.stimInfo.movieNos == movieNumberTimes(1,kk+1));
            if abs(oElecResp.stimInfo.stimAmps(movieInd)) < 2*thresh && abs(oElecResp.stimInfo.stimAmps(nextMovieInd)) >= 2*thresh
                plot([0 100], 0.1*[kk kk], 'k-')
            end
            
            if abs(oElecResp.stimInfo.stimAmps(movieInd)) < 1.5*thresh && abs(oElecResp.stimInfo.stimAmps(nextMovieInd)) >= 1.5*thresh
                plot([0 100], 0.1*[kk kk], 'b-')
            end
        end
        clear movieInd nextMovieInd
    end
    
    disp(['stimulated through ' num2str(max(abs(oElecResp.stimInfo.stimAmps))/thresh) 'x threshold'])
    disp(['(response probability = ' num2str(normcdf(max(abs(oElecResp.stimInfo.stimAmps)), thresh, sig))  ')'])
    
    binEdges = 0:1:100;
    %latenciesMs = spikeTimesCat/20;
    histVals = psthPlotterBase(gca, binEdges, spikeTimesCat);
    title(['n' num2str(special_neuron(ii).id) '(' special_neuron(ii).type ')'])
    set(gca, 'xlim', [-5 100])
    pause
    close(gcf)
end



%% checking to make sure all direct responses of these cells have latencies < 1 ms

neuron(nNormalNeurons+1).elecRespPath = special_neuron(1).originalElecRespPath;
neuron(nNormalNeurons+2).elecRespPath = special_neuron(2).originalElecRespPath;

figure
for ii = 1:length(neuron)
    load(neuron(ii).elecRespPath)
    elecResp = checkForUnfinishedAnalysis(elecResp, 100);
    
    thresh = elecResp.analysis.threshold;
    sig = 1/(elecResp.analysis.erfParams(1)*sqrt(2));

    %amplitude corresponding to histRateMin in cumulative Gaussian fit
    histAmpMin = norminv(0.5, thresh, sig);
    
    aboveThreshMovInd = find(abs(elecResp.stimInfo.stimAmps) > histAmpMin, 1);
    hold on
    for jj = aboveThreshMovInd:length(elecResp.stimInfo.movieNos)
        if ~isempty(elecResp.analysis.type{jj})
             plot_psth(gca, elecResp, elecResp.stimInfo.movieNos(jj))
        end
    end
    hold off
    pause
    cla
end

close(gcf)













