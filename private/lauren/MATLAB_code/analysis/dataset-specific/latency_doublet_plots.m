%plots temporal details about doublets exhibited by ON parasol


cellInfo(1).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data008/';
%cellInfo(1).id = 482;
cellInfo(1).elecRespEarly = 'elecResp_n482_p36_f5.mat';
cellInfo(1).elecRespLate = 'elecResp_n482_p36_f5_ll.mat';
cellInfo(1).visionPath = '2011-01-26-3/data000/data000';
cellInfo(1).id = 482;



binSize = 0.025; %in ms

%seed values for histogram fitting
t0Guess = 0.2;
tauGuess = 0.2;
n = 3;

%checkForUnfinishedAnalysis parameters
nBootstrapReps = 100;
recalcAll = false;

figure
for ii = 1:length(cellInfo);
    
    load([cellInfo(ii).elecRespPath filesep cellInfo(ii).elecRespEarly])
    elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', recalcAll);
    elecResp1 = elecResp; clear elecResp
    
    load([cellInfo(ii).elecRespPath filesep cellInfo(ii).elecRespLate])
    elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', recalcAll);
    elecResp2 = elecResp; clear elecResp
    
    
    nLoneSpikesInt2 = zeros(1,length(elecResp1.stimInfo.movieNos));
    lat2After1Ratio = zeros(1,length(elecResp1.stimInfo.movieNos));
    
    %keeps track of spikes following directly-elicited spikes (spikes
    %within 1ms of stim onset considered directly-elicited)
    lat2After1msNSpikes = zeros(1,length(elecResp1.stimInfo.movieNos));
    lat1msSpikes = zeros(1,length(elecResp1.stimInfo.movieNos));
    
    isiAll = [];
    
    for jj = 1:length(elecResp1.stimInfo.movieNos)
        
        lat1 = elecResp1.analysis.latencies{jj};
        lat2 = elecResp2.analysis.latencies{jj};
        
        %subtract 2 from latencies because stim is applied at sample 2, and convert to ms
        lat1(lat1~=0) = (lat1(lat1~=0)-2)/20;
        lat2(lat2~=0) = (lat2(lat2~=0)-2)/20;
        
        %remove spikes that appear in both intervals (in overlap region) from second interval
        for kk = 1:length(lat2)
            if lat2(kk)~=0 && abs(lat1(kk)-lat2(kk)) < 0.1
                lat2(kk) = 0;
                disp('removed duplicate spike')
            end
        end
        
        nLoneSpikesInt2(jj) = sum(lat2(lat1==0)~=0);
        
        if nLoneSpikesInt2(jj)>0
            disp([ num2str(nLoneSpikesInt2(jj)) ' spikes in second interval without preceding spike in first interval'])
        end
        
        lat2After1msNSpikes(jj) = sum(lat2(lat1~=0 & lat1<1)~=0);
        lat1msSpikes(jj) = sum(lat1~=0 & lat1<1);
        
        lat2After1Ratio(jj) = (sum(lat2(lat1~=0)~=0))/sum(lat2~=0);
        disp([num2str(lat2After1Ratio(jj)*100) ' percent of spikes in second interval follow a spike in the second interval'])
        
        latPairs1 = lat1(lat1~=0 & lat2~=0);
        latPairs2 = lat2(lat1~=0 & lat2~=0);
        
        ISIs = latPairs2 - latPairs1;
        isiAll = [isiAll;ISIs];

        if 0 %look at plots for each stimulus amplitude
            binEdges = 0:binSize:5; %in ms
            hold on
            histVals1 = psthPlotterBase(gca, binEdges, lat1, 'lineColor', [0 0 0], 'fillHist', true);
            histVals2 = psthPlotterBase(gca, binEdges, lat2, 'lineColor', [1 0 0], 'fillHist', true);
            
            %shift time values so that they are in the center of
            %corresponding histogram bins
            
            histVals1(1,:) = histVals1(1,:) + 0.5*binSize;
            histVals2(1,:) = histVals2(1,:) + 0.5*binSize;
            
            %[estParams alpha FWHM ttp fitError] = impRespFitter(histVals1, t0Guess, tauGuess, n, 'makePlot', true, 'plotAxes', gca);
            [estParams alpha FWHM ttp fitError] = impRespFitter(histVals2, t0Guess, tauGuess, n, 'makePlot', true, 'plotAxes', gca);
            
            set(gca, 'xlim', [0 max(binEdges)], 'ylim', [0 length(lat1)/2])
            keyboard
            
            cla
            
            %plot distribution of ISIs
            binEdgesISI = 0:0.1:5;
            psthPlotterBase(gca, binEdgesISI, ISIs, 'lineColor', [0 0 0], 'fillHist', true);
            set(gca, 'xlim', [0 max(binEdgesISI)], 'ylim', [0 10])
            
            keyboard
            cla
        end
    end
    
end

figure
%summary plot of ISIs

binEdgesISI = 0:0.1:5;
isiHistVals = psthPlotterBase(gca, binEdgesISI, isiAll, 'lineColor', [0 0 0], 'fillHist', true);
%shift time values so that they are in the center of
%corresponding histogram bins
isiHistVals(1,:) = isiHistVals(1,:) + 0.5*(isiHistVals(1,2) - isiHistVals(1,1));
[estParams alpha FWHM ttp fitError] = impRespFitter(isiHistVals, t0Guess, tauGuess, n, 'makePlot', true, 'plotAxes', gca);
hold on

%determine peak value of fit
t0 = estParams(1);
tau = estParams(2);
tProj = binEdgesISI(1):0.00001:binEdgesISI(end);
t_sc = (tProj - t0)/tau; %scaled/shifted time values
p = exp(-n*(t_sc-1)).*t_sc.^n;
p(t_sc<0) = 0; %zero out values corresonding to t_xc < 0

fitMax = tProj(p==max(p));

plot(fitMax, max(p)*alpha, 'k*')
text(fitMax+1, max(p)*alpha*1.2, ['peak at t=' num2str(fitMax)])


set(gca, 'xlim', [0 max(binEdgesISI)], 'ylim', [0 20])


%comparison to ACF from WN run
datarun = load_data(cellInfo(1).visionPath);
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
[time acf] = get_correlation(datarun, cellInfo(1).id);
acf = acf(end:-1:1); %reverse along x-axis for plotting



disp([num2str(sum(lat2After1msNSpikes)/sum(lat1msSpikes)*100) '% if spikes occuring within 1ms of stim onset were followed by a spike between ~2 and 4 ms after stim onset'])

