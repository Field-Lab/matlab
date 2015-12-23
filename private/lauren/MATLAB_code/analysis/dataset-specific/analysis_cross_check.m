
clear all

%% cell information
for z=1 %for code folding only
    
    x = 0;
    
    x = x+1;
    cellInfo(x).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002/elecResp_n886_p60';
    cellInfo(x).difficulty = 'easy';
    cellInfo(x).type = 'onMidg';
    cellInfo(x).movieRange = 1:22;
    
    x = x+1;
    cellInfo(x).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/elecResp_n111_p8_w50';
    cellInfo(x).difficulty = 'easy';
    cellInfo(x).type = 'offMidg';
    cellInfo(x).movieRange = 3:2:29;
    
    x = x+1;
    cellInfo(x).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/elecResp_n2_p1_w50';
    cellInfo(x).difficulty = 'easy';
    cellInfo(x).type = 'onPar';
    cellInfo(x).movieRange = 9:2:37;
    
    x = x+1;
    cellInfo(x).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/elecResp_n617_p42_w50';
    cellInfo(x).difficulty = 'easy';
    cellInfo(x).type = 'offPar';
    cellInfo(x).movieRange = 19:2:45;
    
    x = x+1;
    cellInfo(x).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2010-10-28-3/data002/elecResp_n242_p15_w50';
    cellInfo(x).difficulty = 'easy';
    cellInfo(x).type = 'sbc';
    cellInfo(x).movieRange = 41:2:61;
    
    
    x = x+1;
    cellInfo(x).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/elecResp_n302_p21_w50';
    cellInfo(x).difficulty = 'hard';
    cellInfo(x).type = 'onMidg';
    cellInfo(x).movieRange = 29:2:47;
    
    x = x+1;
    cellInfo(x).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/elecResp_n951_p6_w50';
    cellInfo(x).difficulty = 'hard';
    cellInfo(x).type = 'offMidg';
    cellInfo(x).movieRange = 21:2:31;
    
    x = x+1;
    cellInfo(x).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-10-25-4/data002/elecResp_n770_p51_w50';
    cellInfo(x).difficulty = 'hard';
    cellInfo(x).type = 'onPar';
    cellInfo(x).movieRange = 31:2:45;
    
    x = x+1;
    cellInfo(x).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/elecResp_n392_p27_w50';
    cellInfo(x).difficulty = 'hard';
    cellInfo(x).type = 'offPar';
    cellInfo(x).movieRange = 33:2:55;
    
    x = x+1;
    cellInfo(x).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-10-25-4/data002/elecResp_n632_p46_w50';
    cellInfo(x).difficulty = 'hard';
    cellInfo(x).type = 'sbc';
    cellInfo(x).movieRange = 27:2:43;
    
end




%% load analysis

threshAll = zeros(length(cellInfo), 3);
threshAllFraction = zeros(length(cellInfo), 3);

for ii = 1:length(cellInfo)
    
    load([cellInfo(ii).elecRespPath '_lauren']) %analyzed by lauren
    elecRespL = elecResp; clear elecResp
    
    load([cellInfo(ii).elecRespPath '_clare'])
    elecRespC = elecResp; clear elecResp
    
    load([cellInfo(ii).elecRespPath, '_gaby'])
    elecRespG = elecResp; clear elecResp
    
    %make sure data's been fully analyzed in expected movie range only
    for jj = 1:length(elecRespL.stimInfo.movieNos)
        if any(cellInfo(ii).movieRange == elecRespL.stimInfo.movieNos(jj)) %movie should be analyzed
            if isempty(elecRespL.analysis.type{jj}) || ~elecRespL.analysis.finalized(jj)
                disp(['warning: movie ' num2str(elecRespL.stimInfo.movieNos(jj)) ' has not been fully analyzed for ' elecRespL.names.rrs_short_name '_lauren'])
            end
            if isempty(elecRespC.analysis.type{jj}) || ~elecRespC.analysis.finalized(jj)
                disp(['warning: movie ' num2str(elecRespC.stimInfo.movieNos(jj)) ' has not been fully analyzed for ' elecRespC.names.rrs_short_name '_clare'])
            end
            if isempty(elecRespG.analysis.type{jj}) || ~elecRespG.analysis.finalized(jj)
                disp(['warning: movie ' num2str(elecRespG.stimInfo.movieNos(jj)) ' has not been fully analyzed for ' elecRespG.names.rrs_short_name '_gaby'])
            end
        else %movie shouldn't be analyzed
            if ~isempty(elecRespL.analysis.type{jj})
                disp(['warning: movie ' num2str(elecRespL.stimInfo.movieNos(jj)) ' has been analyzed for ' elecRespL.names.rrs_short_name '_lauren'])
            end
            if ~isempty(elecRespC.analysis.type{jj})
                disp(['warning: movie ' num2str(elecRespC.stimInfo.movieNos(jj)) ' has been analyzed for ' elecRespC.names.rrs_short_name '_clare'])
            end
            if ~isempty(elecRespG.analysis.type{jj})
                disp(['warning: movie ' num2str(elecRespG.stimInfo.movieNos(jj)) ' has been analyzed for ' elecRespG.names.rrs_short_name '_gaby'])
            end
        end
    end

    %make sure erf curve fit is current
    elecRespL = checkForUnfinishedAnalysis(elecRespL, 100, 'recalcAll', 0);
    elecResp = elecRespL; save([cellInfo(ii).elecRespPath '_lauren'], 'elecResp'); clear elecResp
    
    elecRespC = checkForUnfinishedAnalysis(elecRespC, 100, 'recalcAll', 0);
    elecResp = elecRespC; save([cellInfo(ii).elecRespPath '_clare'], 'elecResp'); clear elecResp

    elecRespG = checkForUnfinishedAnalysis(elecRespG, 100, 'recalcAll', 0);
    elecResp = elecRespG; save([cellInfo(ii).elecRespPath '_gaby'], 'elecResp'); clear elecResp


    threshAll(ii,1) = elecRespL.analysis.threshold;
    threshAll(ii,2) = elecRespC.analysis.threshold;
    threshAll(ii,3) = elecRespG.analysis.threshold;
    
    threshAllFraction(ii,:) = threshAll(ii,:)/threshAll(ii,1);
    
    clear elecRespL elecRespC elecRespG
end
%%

percentDiffC = 100*abs(threshAllFraction(:,2)-1);
percentDiffG = 100*abs(threshAllFraction(:,3)-1);

disp(['Clare''s thresholds differed from Lauren''s by an average of ' num2str(mean(percentDiffC)),...
    '% (range = ' num2str(min(percentDiffC)) ' - ' num2str(max(percentDiffC))])


disp(['Gaby''s thresholds differed from Lauren''s by an average of ' num2str(mean(percentDiffG))...
    '% (range = ' num2str(min(percentDiffG)) ' - ' num2str(max(percentDiffG))])



















