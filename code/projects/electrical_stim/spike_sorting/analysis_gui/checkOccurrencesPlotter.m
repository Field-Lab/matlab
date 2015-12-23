function checkOccurrencesPlotter(main, elecRespList)

try
    icons = load('icons.mat');
    icons = icons.Icons;
catch
    disp(lasterror)
end

checkOccPlot.figH = figure('position', [500 200 600 800], 'color', 'white', 'Toolbar', 'none',...
    'Menubar', 'none', 'Name', 'check for errors in multiple occurrences of pattern',...
    'KeyPressFcn', @(h_obj, evt)keyCheckOccurrences(h_obj, evt));


%axes
checkOccPlot.axesH = axes('parent', checkOccPlot.figH, 'units', 'pixels',...
    'position', [50 510 220 220]);

checkOccPlot.failAxesH = axes('parent', checkOccPlot.figH, 'units', 'pixels',...
    'position', [290 510 110 220]);

checkOccPlot.eiAxesH = axes('parent', checkOccPlot.figH, 'units', 'pixels',...
    'position', [50 180 300 250]);

%controls
uicontrol(checkOccPlot.figH, 'Style', 'pushbutton', 'String', 'PLOT!',... 
    'Position', [410 700 130 25], 'Callback', @plotNow)


uicontrol(checkOccPlot.figH, 'Style', 'text', 'String', 'amplitude:',...
    'Position', [410 670 100 15], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');

checkOccPlot.ampText = uicontrol(checkOccPlot.figH, 'Style', 'text', 'String', '0',...
    'Position', [410 650 80 15], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');

try
    checkOccPlot.downButton = uicontrol(checkOccPlot.figH,  'Style', 'pushbutton', 'CData', icons.downArrow,...
        'Position', [520 650 24 24], 'Callback', @decreaseAmp);
    checkOccPlot.upButton = uicontrol(checkOccPlot.figH,  'Style', 'pushbutton', 'CData', icons.upArrow,...
        'Position', [550 650 24 24], 'Callback', @increaseAmp);
catch
    checkOccPlot.downButton = uicontrol(checkOccPlot.figH,  'Style', 'pushbutton', 'String', 'down',...
        'Position', [520 652 35 20], 'Callback', @increaseAmp);
    checkOccPlot.upButton = uicontrol(checkOccPlot.figH,  'Style', 'pushbutton', 'String', 'up',...
        'Position', [560 652 35 20], 'Callback', @decreaseAmp);
end


checkOccPlot.preRespCB = uicontrol(checkOccPlot.figH, 'Style', 'checkbox', 'String', 'pre-response?',...
    'position', [410 620 180 25], 'value', false);
checkOccPlot.subYOffsetCB = uicontrol(checkOccPlot.figH, 'Style', 'checkbox', 'String', 'subtract y-offsets? (raw)',...
    'position', [410 595 180 25], 'value', false);

checkOccPlot.plotType = uicontrol(checkOccPlot.figH, 'Style', 'popupmenu', 'String', {'subtract mean','raw data'},...
    'position', [410 560 180 25], 'value', 1, 'BackgroundColor', [1 1 1], 'FontSize', 9);


uicontrol(checkOccPlot.figH, 'Style', 'pushbutton', 'String', 'isolate trace(s)',... 
    'Position', [150 740 100 25], 'Callback', @lassoData)
uicontrol(checkOccPlot.figH, 'Style', 'pushbutton', 'String', 'isolate trace(s)',... 
    'Position', [150 440 100 25], 'Callback', @lassoEIs)


uicontrol(checkOccPlot.figH, 'Style', 'pushbutton', 'String', 'lock displayed analyses',... 
    'Position', [210 60 170 25], 'Callback', @lockAnalyses)

uicontrol(checkOccPlot.figH, 'Style', 'pushbutton', 'String', 'I''m done with this!',... 
    'Position', [210 20 170 25], 'Callback', @cancel)



% elecResp information display
uicontrol(checkOccPlot.figH, 'Style', 'text', 'String', 'displayed elecResp movies:',...
    'Position', [380 460 200 15], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');

checkOccPlot.elecRespText = cell(20, 4);
for ii = 1:20
    checkOccPlot.elecRespText{ii,1} = uicontrol(checkOccPlot.figH, 'Style', 'text', 'String', '',...
        'Position', [380 460-15*ii 160 15], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'fontsize', 9);
    
    checkOccPlot.elecRespText{ii,2} = uicontrol(checkOccPlot.figH, 'Style', 'text', 'String', '',...
        'Position', [540 460-15*ii 18 15], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'fontsize', 9);
    checkOccPlot.elecRespText{ii,3} = uicontrol(checkOccPlot.figH, 'Style', 'text', 'String', '',...
        'Position', [560 460-15*ii 18 15], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'fontsize', 9);
    checkOccPlot.elecRespText{ii,4} = uicontrol(checkOccPlot.figH, 'Style', 'text', 'String', '',...
        'Position', [580 460-15*ii 18 15], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'fontsize', 9);
end


%% 

%store current amplitude

%load elecResp files
elecRespPath = get(main.filePath, 'String');
elecRespsAll = cell(1,length(elecRespList));
for ii = 1:length(elecRespList)
    tmp = load([elecRespPath elecRespList{ii}]);
    elecRespsAll{ii} = tmp.elecResp;
end

%initialize amplitude
currentAmp = elecRespsAll{1}.stimInfo.stimAmps(1);
currentMovies = findCurrMovInds(currentAmp);
set(checkOccPlot.ampText, 'string', currentAmp);
checkOccPlot.plotCurrent = false;

%initialze selected traces
selTraces = [];

%uiwait(checkOccPlot.figH)


    function keyCheckOccurrences(hObject, eventdata) %#ok<INUSL>
        if strcmpi(eventdata.Key, 'uparrow')
            increaseAmp()
            plotNow()
        elseif strcmpi(eventdata.Key, 'downarrow')
            decreaseAmp()
            plotNow()
        elseif strcmpi(eventdata.Key, 'return')
            plotNow()
        elseif strcmpi(eventdata.Key, 'l')
            lockAnalyses()
        end        
    end


    function increaseAmp(~,~)
        if max(currentMovies(1).inds) < length(elecRespsAll{1}.stimInfo.movieNos) %if there are higher movies numbers in first elecResp
            currentAmp = elecRespsAll{1}.stimInfo.stimAmps(max(currentMovies(1).inds)+1);
            currentMovies = findCurrMovInds(currentAmp);
            set(checkOccPlot.ampText, 'string', currentAmp);
        end
        checkOccPlot.plotCurrent = false;
    end

    function decreaseAmp(~,~)
        if min(currentMovies(1).inds) > 1 %if there are higher movies numbers in first elecResp
            currentAmp = elecRespsAll{1}.stimInfo.stimAmps(min(currentMovies(1).inds)-1);
            currentMovies = findCurrMovInds(currentAmp);
            set(checkOccPlot.ampText, 'string', currentAmp);
        end
        checkOccPlot.plotCurrent = false;
    end

    function reloadElecResps()
        %reload elecResps in case they've been updated
        elecRespsAll = cell(1,length(elecRespList));
        for jj = 1:length(elecRespList)
            tmp = load([elecRespPath elecRespList{jj}]);
            elecRespsAll{jj} = tmp.elecResp;
        end
    end

    %check "finalized" status of displayed movies (creates new field for
    %flag)
    function currMovies = checkFinalization(currMovies)
        %checks finalization status of current movies
        for jj = 1:length(elecRespList)
            tmp = load([elecRespPath elecRespList{jj}]);
            currMovies(jj).finalized = [];
            for kk = 1:length(currMovies(jj).inds)
                currMovies(jj).finalized(kk) = tmp.elecResp.analysis.finalized(currMovies(jj).inds(kk));
            end
        end
    end

    function lockAnalyses(~,~)
        if ~checkOccPlot.plotCurrent
            warndlg('plotted data needs to be refreshed')
        else
            for jj = 1:length(elecRespList)
                load([elecRespPath elecRespList{jj}])
                for kk = 1:length(currentMovies(jj).inds)
                    elecResp.analysis.finalized(currentMovies(jj).inds(kk)) = 1; %#ok<STRNU>
                end
                save([elecRespPath elecRespList{jj}], 'elecResp')
            end
            currentMovies = checkFinalization(currentMovies);
            updateElecRespInfoPanel()
        end
    end

    function lassoData(~,~)
        if isfield(checkOccPlot, 'dataTraces')
            selTraces = chooseTracesGui(checkOccPlot.dataTraces);
            plotNow()
        end
    end

    function lassoEIs(~,~)
        if isfield(checkOccPlot, 'eiTraces')
            selEITraces = chooseTracesGui(checkOccPlot.eiTraces);
            successIndeces = find(checkOccPlot.latencies~=0);
            selTraces = successIndeces(selEITraces);
            plotNow()
        end
    end


    function plotNow(hObject, ~)  
        
        %reload elecResps in case they've been updated
        reloadElecResps()
        currentMovies = checkFinalization(currentMovies);
        
        centerChannel = elecRespsAll{1}.cells.recElec;
        dataOrigin = []; %rows correspond with traces in dataTraces, first value specifies elecRespList index, second specifies to movie index
        dataTraces = [];
        latencies = [];
        if exist('hObject', 'var') %if function was called by button press, erase selTraces
            selTraces = [];
        end
        
        plotType = get(checkOccPlot.plotType, 'value');
        
        %load data
        for jj = 1:length(elecRespsAll)            
            for kk = 1:length(currentMovies(jj).inds)
                movieInd = currentMovies(jj).inds(kk);
                if ~isempty(elecRespsAll{jj}.analysis.type{movieInd}) %don't include unanalyzed movies
                    
                    dataTracesTmp=NS_ReadPreprocessedData(elecRespsAll{jj}.names.data_path, '', 0, elecRespsAll{jj}.stimInfo.patternNo,...
                        elecRespsAll{jj}.stimInfo.movieNos(movieInd), 99999);
                    nPulses = size(dataTracesTmp, 1);
                    nSamples = size(dataTracesTmp,3);
                    
                    if nSamples ~= 100 %to line up analysis regions of spatiotemporal probe traces
                        preResp = get(checkOccPlot.preRespCB, 'value');
                        
                        if preResp %analyzes response to prepulse so align traces by first sample
                            startSample = 1;
                        else %analyzes response to stim pulse so align traces by last sample
                            startSample = nSamples-100+1;
                        end
                    else
                        startSample = 1;
                    end
                    
                    dataTracesTmp = squeeze(dataTracesTmp(:,centerChannel,startSample:startSample+59));
                    
                    if plotType == 1 %subtract mean
                        dataTracesTmpMean = mean(dataTracesTmp,1);
                        for i = 1:size(dataTracesTmp,1)
                            dataTracesTmp(i,:) = dataTracesTmp(i,:) - dataTracesTmpMean;
                        end
                    end
                    dataTraces = [dataTraces; dataTracesTmp];
                    dataOriginTmp = [ones(nPulses,1)*jj ones(nPulses,1)*movieInd];
                    dataOrigin = [dataOrigin; dataOriginTmp]; %first index corresponds to elecRespList index, second corresponds to movie index
                    theseLat = elecRespsAll{jj}.analysis.latencies{movieInd};
                    latenciesTmp = zeros(size(theseLat));
                    latenciesTmp(theseLat~=0) = theseLat(theseLat~=0)-startSample+1;
                    latencies = [latencies; latenciesTmp];
                end
                clear dataTracesTmp
            end
        end
                        
        %data manipulation
        subYOffsets = get(checkOccPlot.subYOffsetCB, 'value');
        if subYOffsets
            yOffsets = mean(dataTraces(:,end-35:end),2);
        end
        
%         if plotType == 1 %subtract mean
%             meanTrace = squeeze(mean(dataTraces, 1));
%             if subYOffsets
%                 meanTrace = meanTrace - mean(meanTrace); %set mean value of meanTrace to 0 so that this value isn't subtracted from traces twice
%                 for i = 1:size(dataTraces, 1)
%                     dataTraces(i,:) = dataTraces(i, :) - meanTrace - yOffsets(i); %#ok<AGROW>
%                 end
%             else
%                 for i = 1:size(dataTraces, 1)
%                     dataTraces(i,:) = dataTraces(i, :) - meanTrace; %#ok<AGROW>
%                 end
%             end
%        elseif subYOffsets %raw data with y offsets subtracted
        if plotType ~= 1 && subYOffsets
            for i = 1:size(dataTraces,1)
                dataTraces(i,:) = dataTraces(i, :) - yOffsets(i); %#ok<AGROW>
            end
        end
        
        %plotting (top axes)
        axes(checkOccPlot.axesH) %#ok<MAXES>
        cla; hold on
        plot(dataTraces(latencies==0,:)', 'color', [0.5 0.5 0.5])
        plot(dataTraces(latencies~=0,:)', 'color', [1 0 0])
        
        if ~isempty(selTraces)
            plot(dataTraces(selTraces,:)', 'color', [1 0.5 0], 'LineWidth', 2)
        end
        set(checkOccPlot.axesH, 'xlim', [0 40])
        hold off
        
        axes(checkOccPlot.failAxesH)
        cla; hold on
        plot(dataTraces(latencies==0,:)', 'color', [0 0 0])
        set(checkOccPlot.failAxesH, 'xlim', [5 25], 'ylim', get(checkOccPlot.axesH, 'ylim'))
        hold off
        
        
        %retrieving ei data to plot
        if ~any(latencies > 60)
            centerChannel = elecRespsAll{1}.cells.recElec;
            
            %make a fake elecResp to plug into findStimEi
            elecRespDummy.stimInfo.movieNos = 1;
            elecRespDummy.cells.recElec = 1; %because dataTraces and eiData only contains center electrode's data
            elecRespDummy.cells.goodElecs = 1; %because dataTraces and eiData only contains center electrode's data
            elecRespDummy.cells.main = elecRespsAll{1}.cells.main;
            elecRespDummy.cells.active{1} = [];
            elecRespDummy.stimInfo.nPulses(1) = length(latencies);
            elecRespDummy.analysis.latencies{1} = latencies;
            elecRespDummy.analysis.otherLatencies{1} = [];
            elecRespDummy.analysis.details.analysisFlags{1} = zeros(size(latencies));
                        
            eiData = elecRespsAll{1}.cells.mainEI(centerChannel,:); %ei of target cell
            dataTraces3D = reshape(dataTraces, size(dataTraces,1), 1, size(dataTraces,2));
            
            artFromSubset = get(main.artFromSubsetCheckBox, 'Value');
            %findStimEi(dataTraces, elecResp, movieNumber, tempIndex, eiData, channelsToUse, artFromSubset)
            eiSparse = findStimEi(dataTraces3D, elecRespDummy, 1, 0, {eiData}, 1, artFromSubset);
            targetEiMinPos = find(squeeze(eiData)==min(squeeze(eiData))); %position of minimum on primary electrode
            
            if ~isempty(eiSparse) %if all traces are successes, no ei is calculated
            
                %get rid of zeros in ei
                eiSparse = squeeze(eiSparse);
                eiTraces = eiSparse(latencies~=0,:);
                eiOrigin = dataOrigin(latencies~=0,:);
                
%                 if ~isempty(selTraces)
%                     selTracesBin = false(size(dataTraces,1),1);
%                     selTracesBin(selTraces) = true;
%                     selTracesEI = selTracesBin(latencies~=0);
%                 end
                
                %plotting ei data
                cla(checkOccPlot.eiAxesH)
                axes(checkOccPlot.eiAxesH) %#ok<*MAXES>
                hold on
                
                if subYOffsets
                    plot(eiData(targetEiMinPos-10:targetEiMinPos+15) - mean(eiData(targetEiMinPos-10:targetEiMinPos+15)), 'LineWidth', 2,...
                        'color', 0.5*[1 0 0])
                    for k = 1:size(eiTraces,1)
                        eiTraces(k,:) = eiTraces(k,:)-mean(eiTraces(k,:));
                        plot(eiTraces(k,:), 'color', [1 0 0])
                    end
                else
                    plot(squeeze(eiData(targetEiMinPos-10:targetEiMinPos+15)), 'LineWidth', 2,...
                        'color', 0.5*[1 0 0])
                    for k = 1:size(eiTraces,1)
                        plot(eiTraces(k,:), 'color', [1 0 0])
                    end
                end
                if ~isempty(selTraces)
                    selTracesBin = false(size(dataTraces,1),1);
                    selTracesBin(selTraces) = true;
                    selTracesEI = find(selTracesBin(latencies~=0));
                    for k = 1:length(selTracesEI)
                        plot(eiTraces(selTracesEI(k),:), 'color', [1 0.5 0], 'lineWidth', 2)
                    end
                end
            else
                cla(checkOccPlot.eiAxesH)
                axes(checkOccPlot.eiAxesH)
                text(0.2, 0.5, 'no failure traces.')
                eiTraces = [];
                eiOrigin = [];
            end
        else
            cla(checkOccPlot.eiAxesH)
            axes(checkOccPlot.eiAxesH)
            text(0.2, 0.5, 'Out-of-bound latencies detected.')
            eiTraces = [];
            eiOrigin = [];
        end
                
        %store data to make available to rest of gui
        checkOccPlot.dataTraces = dataTraces;
        checkOccPlot.dataOrigin = dataOrigin;
        checkOccPlot.eiTraces = eiTraces;
        checkOccPlot.eiOrigin = eiOrigin;
        checkOccPlot.latencies = latencies;
        
        updateElecRespInfoPanel()
        checkOccPlot.plotCurrent = true;
    end



    function currMov = findCurrMovInds(targetAmp)
        currMov = struct([]);
        for jj = 1:length(elecRespsAll)
            %find which movie has same amplitude (within 1%)
            matchInd = find(abs((elecRespsAll{jj}.stimInfo.stimAmps - targetAmp)/targetAmp) < 0.005);
            if length(matchInd) < 1
                warndlg([elecRespList{jj} ' doesn''t have at least one matching stimulus amplitude'])
            end
            currMov(jj).inds = matchInd;
        end
        currMov = checkFinalization(currMov);
    end

    function updateElecRespInfoPanel()
        selMovies = checkOccPlot.dataOrigin(selTraces,:);
        
        for jj = 1:20
            for kk = 1:4
                set(checkOccPlot.elecRespText{jj,kk}, 'string', '', 'backgroundcolor', [1 1 1])
            end
        end
        
        for jj = 1:length(elecRespList)
            undInd = strfind(elecRespList{jj}, '_');
            matInd = strfind(elecRespList{jj}, '.mat');
            set(checkOccPlot.elecRespText{jj,1}, 'string', elecRespList{jj}(undInd(1)+1:matInd-1))
            for kk = 1:length(currentMovies(jj).inds)
                movieInd = currentMovies(jj).inds(kk);
                if any((selMovies(:,1) == jj).*(selMovies(:,2)==movieInd)) %at least one selected trace from this movie
                    set(checkOccPlot.elecRespText{jj,kk+1}, 'string', num2str(movieInd),...
                        'ForegroundColor', [1 0.5 0], 'FontWeight', 'bold')
                    if isempty(elecRespsAll{jj}.analysis.type{movieInd})
                        error('selected trace appears to be from an unanalyzed movie - error in code')
                    end
                elseif ~isempty(elecRespsAll{jj}.analysis.type{movieInd}) %only display for analyzed movies
                    set(checkOccPlot.elecRespText{jj,kk+1}, 'string', num2str(movieInd),...
                        'ForegroundColor', [0 0 0], 'FontWeight', 'normal')
                end
                
                %display as bold if movies have been finalized
                if currentMovies(jj).finalized(kk)
                    set(checkOccPlot.elecRespText{jj,kk+1}, 'backgroundcolor', [1 1 0.5])
                else
                    set(checkOccPlot.elecRespText{jj,kk+1}, 'backgroundcolor', [1 1 1])
                end
            end
        end
    end



    function cancel(~, ~)
        close(checkOccPlot.figH)
    end
end