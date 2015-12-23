function analysisGUI()

%June 2011: mades some changes to how stimulus is plotted (within
%funH.showStimulus) to deal with situation when number of electrodes
%changes from movie to movie - if number of electrodes in current stimulus
%is not consistent with elecResp.stimInfo.electrodes, full stimulus is 
%plotted but electrode identities are not displayed

clear all
close all

elecResp = [];

funH.switchMode = @switchMode;
funH.resetForNewData = @resetForNewData;
funH.decreaseMovieNo = @decreaseMovieNo;
funH.increaseMovieNo = @increaseMovieNo;
funH.refreshAllDummy = @refreshAllDummy;
funH.reanalyze = @reanalyze;
funH.tempAnalyze = @tempAnalyze;
funH.showStimulus = @showStimulus;
funH.fitToErf = @fitToErf;
funH.findNeuronsOnElecs = @findNeuronsOnElecs;
funH.addActiveNeuronDummy = @addActiveNeuronDummy;
funH.createElecRespDummy = @createElecRespDummy;
funH.removeFalsePos = @removeFalsePos;
funH.addMissedSpike = @addMissedSpike;
funH.clearAnalysis = @clearAnalysis;
funH.clearRemainingAnalysis = @clearRemainingAnalysis;
funH.clearPrecedingAnalysis = @clearPrecedingAnalysis;
funH.lockAnalysis = @lockAnalysis;
funH.refreshEiPlotsDummy = @refreshEiPlotsDummy;
funH.expandAxes = @expandAxes;
funH.otherElecsSwitch = @otherElecsSwitch;
funH.analyzeMultiple = @analyzeMultipleOcc;
funH.switchLimMode = @switchLimMode;
funH.altPlots = @altPlots;
funH.displayOtherTemplates = @dispOtherTemplates;
funH.increaseOccNo = @increaseOccNo;

main = generateAnalysisGuiLayout(funH);

%% initialization tasks

set(main.gui, 'Visible', 'on')

elecResp.dummyField = [];

set(main.gui, 'KeyPressFcn', @(h_obj, evt)keyCheck(h_obj, evt));

%% callbacks

    function keyCheck(hObject, eventdata) %#ok<INUSL>
        %disp(eventdata.Key)
        %class(eventdata.Key)
        if strcmpi(eventdata.Key, 'uparrow')
            increaseMovieNo()
            main.movieDirection = 'up';
            [main, elecResp] = refreshAll(main);
        elseif strcmpi(eventdata.Key, 'downarrow')
            decreaseMovieNo()
            main.movieDirection = 'down';
            [main, elecResp] = refreshAll(main);
        elseif strcmpi(eventdata.Key, 'return')
            [main, elecResp] = refreshAll(main);
        elseif strcmpi(eventdata.Key, 'l')
            lockAnalysis()
        elseif strcmpi(eventdata.Key, 'g')
            reanalyze()
        elseif strcmpi(eventdata.Key, 'r')
            removeFalsePos()
        end
    end

    function refreshAllDummy(hObject, eventdata) %#ok<INUSD>
        main.movieDirection = '';
        [main, elecResp] = refreshAll(main);
    end

    function refreshEiPlotsDummy(hObject, eventdata) %#ok<INUSD>
        currentMovie = str2double(get(main.movieNo, 'String'));
        if ~sum(elecResp.stimInfo.movieNos == currentMovie)
            warnH = warndlg('Current movie number is not valid for this dataset/pattern.');
            uiwait(warnH)
            return
        end

        if exist([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo) filesep 'p' num2str(elecResp.stimInfo.patternNo) '_m' num2str(currentMovie)], 'file')
            dataTraces=NS_ReadPreprocessedData([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo)],...
                '', 0, elecResp.stimInfo.patternNo, currentMovie, 99999);
            refreshEiPlots(main, elecResp, dataTraces)
        elseif exist([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo) '_m' num2str(currentMovie)], 'file')
            dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo,...
                currentMovie, 99999);
            refreshEiPlots(main, elecResp, dataTraces)
        else
            warnH = warndlg([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo) '_m' num2str(currentMovie) ' couldn''t be found']);
            uiwait(warnH)
        end
    end

    function increaseMovieNo(hObject, eventdata) %#ok<INUSD>
        currentMovie = str2double(get(main.movieNo, 'String'));
        if isfield(elecResp, 'stimInfo')
            movies = elecResp.stimInfo.movieNos;
            movieIndex = find(movies == currentMovie, 1);
            if movieIndex < length(movies)
                set(main.movieNo, 'String', movies(movieIndex+1))
            elseif isempty(movieIndex)
                warnH = warndlg('Current movie number is not valid for this dataset/pattern. Resetting to first movie number...');
                uiwait(warnH)
                set(main.movieNo, 'String', movies(1))
            end
        else
            set(main.movieNo, 'String', currentMovie + 1)
        end
        main.tempAnalysisNoSave = [];
    end

    function decreaseMovieNo(hObject, eventdata) %#ok<INUSD>
        currentMovie = str2double(get(main.movieNo, 'String'));
        if isfield(elecResp, 'stimInfo')
            movies = elecResp.stimInfo.movieNos;
            movieIndex = find(movies == currentMovie, 1);
            if movieIndex > 1
                set(main.movieNo, 'String', movies(movieIndex-1))
            elseif isempty(movieIndex)
                warnH = warndlg('Current movie number is not valid for this dataset/pattern. Resetting to first movie number...');
                uiwait(warnH)
                set(main.movieNo, 'String', movies(1))
            end
        elseif currentMovie > 1
            set(main.movieNo, 'String', currentMovie - 1)
        end
        main.tempAnalysisNoSave = [];
    end

    function increaseOccNo(~, ~)
        erName = get(main.fileName, 'String');
        o_ind = strfind(erName, 'o');
        if length(o_ind) == 1
            after_o_ind = [strfind(erName, '.') strfind(erName, '_')];
            after_o_ind = min(after_o_ind(after_o_ind>o_ind));
            if ~isempty(after_o_ind)
                o_val = str2double(erName(o_ind+1:after_o_ind-1));
                set(main.fileName, 'String', [erName(1:o_ind) num2str(o_val+1) erName(after_o_ind:end)])
            else
                o_val = str2double(erName(o_ind+1:end));
                set(main.fileName, 'String', [erName(1:o_ind) num2str(o_val+1)])
            end
        elseif isempty(o_ind)
            warndlg('no ''o'' in elecResp name')
        else
            warndlg('more than one ''o'' in elecResp name -- if you want to use this button, make the code smarter!!')
        end
    end

    function analyzeMultipleOcc(hObject, eventdata) %#ok<INUSD>
        checkMultipleOccurrences(main)
    end

    function reanalyze(hObject, eventdata) %#ok<INUSD>
        
        filePath = get(main.filePath, 'String');
        fileName = get(main.fileName, 'String');
        
        temp = load([filePath fileName]);
        elecResp = temp.elecResp;
        
%         if ~isfield(elecResp, 'stimInfo')
%             warnH = warndlg('need to load up an elecResp file first!');
%             uiwait(warnH)
%             return
%         end

        %warns user if analysis was previously finalized (set as 'locked')
        currentMovie = str2double(get(main.movieNo, 'String'));
        movieIndex = find(elecResp.stimInfo.movieNos == currentMovie);
        if elecResp.analysis.finalized(elecResp.stimInfo.movieNos == currentMovie)
            response = lockOverwriteGui();
            if ~response
                return
            end
        end
        
        modelTypeValue = get(main.reanalyzeType, 'Value');
        switch modelTypeValue
            case 1
                modelType = 'linkage';
            case 2
                modelType = 'prevArtifact';
            case 3
                modelType = 'nextArtifact';
            case 4
                modelType = 'ttx';
            case 5
                modelType = 'handPicked';
            case 6
                modelType = 'currentArtifact';
            case 7
                modelType = 'otherOccurrences';
        end
        
        params.saveFigures = 0;
        
        if get(main.fixedLimsCB, 'value') %use limits defined in elecResp file
            if ~isempty(elecResp.analysis.details.residCalcWindow{movieIndex})
                
                spikeMinWindow = convertTempOffsetWindow(temp.elecResp, movieIndex);

                %updates gui values
                %set(main.xAxisLimitLow, 'string', num2str(elecResp.analysis.details.residCalcWindow{movieIndex}(1)-10));
                %set(main.xAxisLimitHigh, 'string', num2str(elecResp.analysis.details.residCalcWindow{movieIndex}(2)+5))
                set(main.residualLimitLow, 'string', num2str(elecResp.analysis.details.residCalcWindow{movieIndex}(1)))
                set(main.residualLimitHigh, 'string', num2str(elecResp.analysis.details.residCalcWindow{movieIndex}(2)))
                set(main.offsetLimitLow, 'string', num2str(spikeMinWindow(1)))
                set(main.offsetLimitHigh, 'string', num2str(spikeMinWindow(2)))
                set(main.offsetStep, 'string', num2str(elecResp.analysis.details.tempOffsetStep{movieIndex}))
            end
        end
        
        params.shiftLimSpikeMin = [str2double(get(main.offsetLimitLow, 'String'))...
            str2double(get(main.offsetLimitHigh, 'String'))];
        params.shiftStep = str2double(get(main.offsetStep, 'String'));
        params.residLim = [str2double(get(main.residualLimitLow, 'String'))...
            str2double(get(main.residualLimitHigh, 'String'))];
        params.spikeMinExclusionReg = str2double(get(main.spikeMinExclusionReg, 'String'));

        
        analyzeMultiple = get(main.multipleMoviesCheckBox, 'Value');
        clearFlags = get(main.clearFlagsCheckBox, 'Value');
        noRefit = get(main.noRefitCheckBox, 'Value');
        
        if analyzeMultiple
            movieNos = chooseMoviesGui(elecResp);
        else
            movieNos = currentMovie;
        end
        
        if clearFlags
            analysisFlags = zeros(elecResp.stimInfo.nPulses(movieIndex),...
                length(elecResp.cells.active{movieIndex})+1);
        end
        
        
        for i = 1:length(movieNos)
            if ~clearFlags
                analysisFlags = elecResp.analysis.details.analysisFlags{elecResp.stimInfo.movieNos == movieNos(i)};
            end
            
            elecResp = templateMatchClustering(elecResp, movieNos(i), params, 'modelType', modelType,...
                'progressBar', 'off', 'analysisFlags', analysisFlags, 'noRefit', logical(noRefit));

            filePath = get(main.filePath, 'String');
            fileName = get(main.fileName, 'String');
            save([filePath filesep fileName], 'elecResp')
            
            set(main.movieNo, 'String', num2str(movieNos(i)))

            [main, elecResp] = refreshAll(main);

            disp(['done analyzing movie ' num2str(movieNos(i))])
        end
    end

    function addActiveNeuronDummy(hObject, eventdata) %#ok<INUSD>
        if get(main.analyzeModeButton, 'Value')
        
            if ~isfield(elecResp, 'stimInfo')
                warnH = warndlg('need to load up an elecResp file first!');
                uiwait(warnH)
            else
                currentMovie = str2double(get(main.movieNo, 'String'));
                movieIndex = find(elecResp.stimInfo.movieNos == currentMovie, 1);

                [elecResp action] = addActiveNeuron(elecResp, movieIndex);

                if ~strcmpi(action, 'cancel')
                    filePath = get(main.filePath, 'String');
                    fileName = get(main.fileName, 'String');
                    save([filePath filesep fileName], 'elecResp')

                    reanalyze()
                end
            end
        else %explore mode: create new external ei file
            elecRespDummy.cells.recElec = str2double(get(main.centerElec, 'String'));
            elecRespDummy.stimInfo.movieNos = str2double(get(main.movieNo, 'String'));
            elecRespDummy.names.data_path = get(main.filePath, 'String');
            elecRespDummy.stimInfo.patternNo = str2double(get(main.patternNo, 'String'));        
            patternNo = str2double(get(main.patternNo, 'String'));
            movieNo = str2double(get(main.movieNo, 'String'));
            
            if isnan(elecRespDummy.stimInfo.patternNo) && ~isempty(main.tempAnalysisNoSave)
                elecRespDummy.analysis.latencies{1} = main.tempAnalysisNoSave;
                elecRespDummy.stimInfo.patternNo = get(main.patternNo, 'String');
            elseif size(main.tempAnalysis, 1) >= patternNo && size(main.tempAnalysis, 2) >= movieNo && ~isempty(main.tempAnalysis{patternNo, movieNo})
                elecRespDummy.analysis.latencies{1} = main.tempAnalysis{patternNo, movieNo};
            else
                error('need to (temporarily) analyze data before creating ei')
            end
            
            xLimits = [str2double(get(main.xAxisLimitLow, 'String')) str2double(get(main.xAxisLimitHigh, 'String'))];

            elecRespDummy.stimInfo.nPulses(1) = length(elecRespDummy.analysis.latencies{1});
            
            ei = generateEIFromStimData(elecRespDummy, 1, 'xLimits', xLimits); %movieIndex = 1 because dummy elecResp only has one movie
            
            paramsPath = get(main.eiPath, 'String');
            paramsPath = [paramsPath(1:end-2) 'params'];

            [newID action] = chooseNeuronID(paramsPath);
            
            if ~isempty(ei) && strcmp(action, 'success'); %new ei successfully generated
                uisave('ei', ['eiFile' num2str(newID), '.mat']);
            end
            
        end
    end

    function lockAnalysis(hObject, eventdata) %#ok<INUSD>
        if ~isfield(elecResp, 'stimInfo')
            warnH = warndlg('need to load up an elecResp file first!');
            uiwait(warnH)
        else
            filePath = get(main.filePath, 'String');
            fileName = get(main.fileName, 'String');
            
            %reloads elecResp file to avoid problems with interrupting 'refreshAll'
            temp = load([filePath fileName]);
            elecResp = temp.elecResp;
            
            currentMovie = str2double(get(main.movieNo, 'String'));
            movieIndex = find(elecResp.stimInfo.movieNos == currentMovie, 1);
            elecResp.analysis.finalized(movieIndex) = 1;
            set(main.lockLatenciesButton, 'enable', 'off', 'String', 'analysis locked')

            save([filePath filesep fileName], 'elecResp')
        end
    end

    function clearAnalysis(hObject, eventdata, movieIndexArg) %#ok<INUSL>ex
        filePath = get(main.filePath, 'String');
        fileName = get(main.fileName, 'String');
        
        temp = load([filePath fileName]);
        elecResp = temp.elecResp;
        
        if exist('movieIndexArg', 'var')
            movieIndex = movieIndexArg;
        else
            currentMovie = str2double(get(main.movieNo, 'String'));
            movieIndex = find(elecResp.stimInfo.movieNos == currentMovie);
        end
        
        if elecResp.analysis.finalized(movieIndex)
            response = lockOverwriteGui();
            if ~response
                return
            end
        end
        
        
        elecResp.analysis.type{movieIndex}          = [];
        elecResp.analysis.estArtifact{movieIndex}   = [];
        elecResp.analysis.latencies{movieIndex}     = [];
        elecResp.analysis.otherLatencies{movieIndex} = [];
        elecResp.analysis.finalized(movieIndex)      = 0;
        elecResp.analysis.successRates(movieIndex)   = 0;
        elecResp.analysis.erfCurrent                   = 0;
        
        elecResp.analysis.details.residCalcWindow{movieIndex}  = [];
        elecResp.analysis.details.tempOffsetWindow{movieIndex} = [];
        elecResp.analysis.details.tempOffsetStep{movieIndex}   = [];
        elecResp.analysis.details.linkCalcWindow{movieIndex}   = [];
        elecResp.analysis.details.linkThresh{movieIndex}       = [];
        elecResp.analysis.details.nBranches{movieIndex}        = [];
        elecResp.analysis.details.analysisFlags{movieIndex}    =...
            zeros(elecResp.stimInfo.nPulses(movieIndex), length(elecResp.cells.active{movieIndex})+1);
        
        filePath = get(main.filePath, 'String');
        fileName = get(main.fileName, 'String');
        save([filePath filesep fileName], 'elecResp')
        [main, elecResp] = refreshAll(main);
    end

    function clearRemainingAnalysis(hObject, eventdata) %#ok<INUSD>
        currentMovie = str2double(get(main.movieNo, 'String'));
        movieNos = elecResp.stimInfo.movieNos;
        movieIndex = find(movieNos == currentMovie);
        
        for i = movieIndex:length(movieNos)
            clearAnalysis('','', i)
        end
    end

    function clearPrecedingAnalysis(hObject, eventdata) %#ok<INUSD>
        currentMovie = str2double(get(main.movieNo, 'String'));
        movieNos = elecResp.stimInfo.movieNos;
        movieIndex = find(movieNos == currentMovie);
        
        for i = 1:movieIndex
            clearAnalysis('','', i)
        end
    end


    function showStimulus(hObject, eventdata) %#ok<INUSD>
        
        currentMovie = str2double(get(main.movieNo, 'String'));
        
        if get(main.analyzeModeButton, 'Value')        
            if ~isfield(elecResp, 'stimInfo')
                warnH = warndlg('need to load up an elecResp file first!');
                uiwait(warnH)
            else
                movieIndex = find(elecResp.stimInfo.movieNos == currentMovie);
                if ~isempty(movieIndex)
                    nElec = length(elecResp.stimInfo.electrodes);
                    nElecActual = size(elecResp.stimInfo.pulseVectors{movieIndex},1); %because nElec might not reflect actual number of electrodes involved in cases where this differs from movie to movie
                    stimFig = figure('Position', [100 500 (200*nElecActual+60) 230], 'Color', 'white');
                    maxAmp = max(max(abs(elecResp.stimInfo.pulseVectors{movieIndex}(:,2,:))));
                    for i = 1:nElecActual
                        axes('Parent', stimFig, 'Units', 'pixels', 'Position', [(60+200*(i-1)) 40 150 150])
                        plot(squeeze(elecResp.stimInfo.pulseVectors{movieIndex}(i,1,:)),...
                            squeeze(elecResp.stimInfo.pulseVectors{movieIndex}(i,2,:)));
                        if nElec == nElecActual
                            title(['electrode ' num2str(elecResp.stimInfo.electrodes(i))], 'FontSize', 14, 'FontWeight', 'bold')
                        end
                        xlabel('time (\mus)')
                        set(gca, 'ylim', [-maxAmp*1.1 maxAmp*1.1])
                        if i == 1
                            ylabel('current (\muA)')
                        end
                    end
                else
                    warnH = warndlg('Movie number is not valid.');
                    uiwait(warnH)
                end
            end
        else %explore mode
            
            dataPath = get(main.filePath, 'String');
            patternNo = get(main.patternNo, 'String');
            interval = 50; %length of time before stimulus is updated, in microseconds

      
            
            try
                [amps electrodes stimAmpVectors] = getStimAmps(dataPath, patternNo, currentMovie);
            catch
                disp(lasterror)
                disp('unable to get stimulus information')
                return
            end
            
            plotMax = max(max(stimAmpVectors)) + 0.02;
            plotMin = min(min(stimAmpVectors)) - 0.02;
            
            nElec = length(electrodes);
            patternLength = size(stimAmpVectors, 2);
            
            stimFig = figure('Position', [100 500 (200*nElec+60) 230], 'Color', 'white');
            for j = 1:nElec
                pulseVectors = zeros(nElec, 2, (patternLength-1)*2);
                for k = 1:patternLength-1
                    pulseVectors(j, 1, 2*(k-1)+2) = k*interval;
                    pulseVectors(j, 1, 2*(k-1)+3) = k*interval;
                    pulseVectors(j, 2, 2*(k-1)+1) = stimAmpVectors(j,k);
                    pulseVectors(j, 2, 2*(k-1)+2) = stimAmpVectors(j,k);
                end
                pulseVectors(j, 2, 2*(patternLength-1)+1) = stimAmpVectors(j,patternLength);
                
                axes('Parent', stimFig, 'Units', 'pixels', 'Position', [(60+200*(j-1)) 40 150 150])
                plot(squeeze(pulseVectors(j,1,:)), squeeze(pulseVectors(j,2,:)));
                title(['electrode ' num2str(electrodes(j))], 'FontSize', 14, 'FontWeight', 'bold')
                xlabel('time (\mus)')
                set(gca, 'ylim', [plotMin plotMax])
                if j == 1
                    ylabel('current (\muA)')
                end
            end
        end
    end


    function createElecRespDummy(hObject, eventdata) %#ok<INUSD>
        createElecRespStructGui(main)
    end


    function findNeuronsOnElecs(hObject, eventdata) %#ok<INUSD> 
        if get(main.analyzeModeButton, 'Value')
            if ~isfield(elecResp, 'stimInfo')
                warnH = warndlg('need to load up an elecResp file first! better luck next time.');
                uiwait(warnH)
            elseif get(main.exploreOtherElecsCheckBox, 'value')
                if strcmp('/', filesep)
                    elecRespFake.cells.recElec = str2double(get(main.centerExploreNoSpikes, 'String'));
                    %elecRespFake.cells.goodElecs = elecResp.cells.goodElecs;
                    elecRespFake.names.rrs_ei_path = elecResp.names.rrs_ei_path;
                    elecRespFake.names.rrs_params_path = elecResp.names.rrs_params_path;
                    findNeuronsOnElecsGui(elecRespFake, main)
                end
            else
                findNeuronsOnElecsGui(elecResp, main)
            end
        else
            elecRespFake.names.rrs_ei_path = get(main.eiPath, 'String');
            if exist(elecRespFake.names.rrs_ei_path, 'file') && strcmp(elecRespFake.names.rrs_ei_path(end-2:end), '.ei')
                elecRespFake.cells.recElec = str2double(get(main.centerElec, 'String'));
                %elecRespFake.cells.goodElecs = [];
                elecRespFake.names.rrs_params_path = [elecRespFake.names.rrs_ei_path(1:end-2) 'params'];
                findNeuronsOnElecsGui(elecRespFake, main)
            else
                warnH = warndlg('specified ei file path is invalid.  get your act together!');
                uiwait(warnH)
            end
        end
    end


    function expandAxes(hObject, eventdata) %#ok<INUSD>
        currentMovie = str2double(get(main.movieNo, 'String'));
        
        if get(main.analyzeModeButton, 'Value')
            dataPath = elecResp.names.data_path;
            pattern = elecResp.stimInfo.patternNo;
            if ~(elecResp.stimInfo.movieNos == currentMovie)
                warnH = warndlg('movie number is not valid for this elecResp file.');
                uiwait(warnH)
                return
            end

        else
            dataPath = get(main.filePath, 'String');
            pattern = get(main.patternNo, 'String');
        end
        
        if exist([dataPath filesep 'p' pattern filesep 'p' pattern '_m' num2str(currentMovie)], 'file')
            dataTraces=NS_ReadPreprocessedData([dataPath filesep 'p' pattern], '', 0, pattern, currentMovie, 99999);
        elseif exist([dataPath filesep 'p' pattern '_m' num2str(currentMovie)], 'file')
            dataTraces=NS_ReadPreprocessedData(dataPath, '', 0, pattern, currentMovie, 99999);
        else
            warnH = warndlg([dataPath filesep 'p' pattern '_m' num2str(currentMovie) ' couldn''t be found']);
            uiwait(warnH)
            return
        end
        
        if strcmp('/', filesep)
            for i = 1:7
                for j = 1:2
                    if hObject == main.aButton{j,i}
                        plotElecIndeces = [j i];
                    end
                end
            end
        else
            for i = 1:7
                if hObject == main.aButton{1,i}
                    plotElecIndeces = [1 i];
                end
            end
        end
        
        if strcmp('/', filesep)
            newFig = figure('position', [400 400 600 400]');
        else
            newFig = figure('position', [200 200 600 400]');
        end
        plotAxes = axes('parent', newFig);
        refreshEiPlots(main, elecResp, dataTraces, 'elecOverride', plotElecIndeces, 'axesOverride', plotAxes)
    end

    function dispOtherTemplates(~,~)
        tmp = dispOtherTemplatesGui(main);
        if isempty(tmp) || tmp ~= 0
            main.otherDispEIs = tmp;
        end
        main = refreshAll(main);
    end


    function addMissedSpike(hObject, eventdata) %#ok<INUSD>
        
        filePath = get(main.filePath, 'String');
        fileName = get(main.fileName, 'String');
        temp = load([filePath fileName]);
        elecResp = temp.elecResp;
        
        movies = elecResp.stimInfo.movieNos;
        currentMovie = str2double(get(main.movieNo, 'String'));
        movieIndex = find(movies == currentMovie);
        latencies = elecResp.analysis.latencies{movieIndex};
        centerChannel = elecResp.cells.recElec;
        plotType = get(main.dataDisplayType, 'Value');
        
        if elecResp.analysis.finalized(movieIndex)
            response = lockOverwriteGui();
            if ~response
                return
            end
        end
        
        
        if exist([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo) filesep 'p' num2str(elecResp.stimInfo.patternNo) '_m' num2str(currentMovie)], 'file')
            dataTraces=NS_ReadPreprocessedData([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo)],...
                '', 0, elecResp.stimInfo.patternNo, currentMovie, 99999);
        elseif exist([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo) '_m' num2str(currentMovie)], 'file')
            dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo,...
                currentMovie, 99999);
        else
            warnH = warndlg([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo) '_m' num2str(currentMovie) ' couldn''t be found']);
            uiwait(warnH)
            return
        end
        
        if plotType == 2 %subtract mean
            subtractionVector = squeeze(mean(dataTraces(:, centerChannel, :), 1));
        elseif plotType == 3 %subtract estimated artifact
            if isempty(elecResp.analysis.estArtifact{movieIndex})
                warnh = warndlg('No artifact estimate has been calculated yet.  Displaying raw data.');
                uiwait(warnh)
                subtractionVector = zeros(size(dataTraces, 3), 1);
            else
                subtractionVector = elecResp.analysis.estArtifact{movieIndex}(elecResp.cells.goodElecs == centerChannel, :)';
            end
        end

        dataToPlot = zeros(size(dataTraces, 1), size(dataTraces, 3));
        if any(plotType == [2 3]) %subtract mean or estimated artifact
            for i = 1:elecResp.stimInfo.nPulses(movieIndex)
                dataToPlot(i, :) = squeeze(dataTraces(i, centerChannel, :)) - subtractionVector;
            end
        else %display raw data
            dataToPlot = squeeze(dataTraces(:, centerChannel, :));
        end
        
        if ~isempty(elecResp.analysis.otherLatencies{movieIndex})
            initialSuccesses = sum(elecResp.analysis.otherLatencies{movieIndex}, 2) + latencies; %vector where 0 means none of the neurons spiked
        else
            initialSuccesses = latencies;
        end
        
        if all(initialSuccesses)
            warnh = warndlg('No failure traces.  If spike is recognized as coming from wrong neuron, use "remove false positive" and reanalyze for other neuron templates.');
            uiwait(warnh)
            return
        end
                
        if size(dataToPlot,2) > max([80 diff(elecResp.analysis.details.residCalcWindow{movieIndex})+30]) %trace is too long to display well
            xLim = elecResp.analysis.details.residCalcWindow{movieIndex} + [-20 20];
        else
            xLim = [];
        end
        
        traceIndeces = chooseTracesGui(dataToPlot, 'initialSelectedBin', ~initialSuccesses, 'returnSingleIndex', 'false', 'xPlotLim', xLim);
        
        for i = 1:length(traceIndeces)
            elecResp = templateMatchClustering1Trace(elecResp, currentMovie, traceIndeces(i),...
                3*ones(1, length(elecResp.cells.active{movieIndex})+1));
        end
        
        filePath = get(main.filePath, 'String');
        fileName = get(main.fileName, 'String');
        save([filePath filesep fileName], 'elecResp')

        [main, elecResp] = refreshAll(main);
        disp('done analyzing')
    end

    function removeFalsePos(hObject, eventdata) %#ok<INUSD>
        filePath = get(main.filePath, 'String');
        fileName = get(main.fileName, 'String');
        temp = load([filePath fileName]);
        elecResp = temp.elecResp;
        
        movies = elecResp.stimInfo.movieNos;
        currentMovie = str2double(get(main.movieNo, 'String'));
        movieIndex = find(movies == currentMovie);
        nTemplates = length(elecResp.cells.active{movieIndex}) + 1;
        
        if elecResp.analysis.finalized(movieIndex)
            response = lockOverwriteGui();
            if ~response
                return
            end
        end
        
        [template, channelToPlot, reanalyze] = chooseTemplateAndChannelGui(elecResp, movieIndex);
        
        if template == 0 %user hit cancel button
            return
        elseif template == elecResp.cells.main
            tempIndex = 0; %signifies that main neuron was selected
        else
            tempIndex = find(elecResp.cells.active{movieIndex} == template);
        end
        
        
        
        % load eiData
        %eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(elecResp.names.rrs_ei_path);
        eiData = cell(nTemplates, 1); %stores EIs of active neurons: eiData{neuronIndex}(channels, samples)
        %eiData{1} = eiFile.getImage(elecResp.cells.main);
        %eiData{1} = squeeze(eiData{1}(1, 2:end, :));
        eiData{1} = elecResp.cells.mainEI;
        for i = 2:nTemplates
            %eiData{i} = eiFile.getImage(elecResp.cells.active{movieIndex}(i-1));
            %eiData{i} = squeeze(eiData{i}(1, 2:end, :));
            eiData{i} = elecResp.cells.allEIs{elecResp.cells.all == elecResp.cells.active{movieIndex}(i-1)};
        end
        clear eiFile
        
        %load dataTraces
        if exist([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo) filesep 'p' num2str(elecResp.stimInfo.patternNo) '_m' num2str(currentMovie)], 'file')
            dataTraces=NS_ReadPreprocessedData([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo)],...
                '', 0, elecResp.stimInfo.patternNo, currentMovie, 99999);
        elseif exist([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo) '_m' num2str(currentMovie)], 'file')
            dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo,...
                currentMovie, 99999);
        else
            warnH = warndlg([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo)...
                '_m' num2str(currentMovie) ' couldn''t be found']);
            uiwait(warnH)
            return
        end

        if tempIndex == 0
            spikeBin = elecResp.analysis.latencies{movieIndex}~=0;
        else
            spikeBin = elecResp.analysis.otherLatencies{movieIndex}(:,tempIndex)~=0;
        end
        
        if ~sum(spikeBin)
            warnh = warndlg('No spikes to remove for this neuron.');
            uiwait(warnh)
            return
        end
        
        artFromSubset = get(main.artFromSubsetCheckBox, 'Value');
        dataToPlot = findStimEi(dataTraces, elecResp, currentMovie, tempIndex, eiData, channelToPlot, artFromSubset);
        dataToPlot = reshape(dataToPlot, size(dataToPlot, 1), []);
        
        traceIndeces = chooseTracesGui(dataToPlot, 'returnSingleIndex', 'false', 'initialSelectedBin', spikeBin);
        
        if reanalyze
            for i = 1:length(traceIndeces)
                flags = zeros(1, nTemplates);
                if template == elecResp.cells.main
                    flags(1) = 4;
                else
                    flags(find(elecResp.cells.active{movieIndex}==template)+1) = 4;
                end
                elecResp = templateMatchClustering1Trace(elecResp, currentMovie, traceIndeces(i), flags);
            end
        else
            for i = 1:length(traceIndeces)
                elecResp.analysis.details.analysisFlags{movieIndex}(traceIndeces(i),:) = ones(1, nTemplates);
                elecResp.analysis.details.analysisFlags{movieIndex}(traceIndeces(i), tempIndex+1) = 4;
                if tempIndex == 0
                    elecResp.analysis.latencies{movieIndex}(traceIndeces(i)) = 0;
                else
                    elecResp.analysis.otherLatencies{movieIndex}(traceIndeces(i), tempIndex) = 0;
                end
                elecResp.analysis.finalized(movieIndex) = 0;
            end
            
            if tempIndex == 0
                elecResp.analysis.successRates(movieIndex) = ...
                    sum(elecResp.analysis.latencies{movieIndex}~=0)/elecResp.stimInfo.nPulses(movieIndex);
            end
        end
        filePath = get(main.filePath, 'String');
        fileName = get(main.fileName, 'String');
        save([filePath '/' fileName], 'elecResp')
        [main, elecResp] = refreshAll(main);
        reanalyze();
    end

    function tempAnalyze(hObject, eventdata) %#ok<INUSD> %hand clustering in explore mode
        dataPath = get(main.filePath, 'String');
        patternNo = str2double(get(main.patternNo, 'String'));
        movieNo = str2double(get(main.movieNo, 'String'));
        centerChannel = str2double(get(main.centerElec, 'String'));
        dataTraces=NS_ReadPreprocessedData(dataPath, '', 0, get(main.patternNo, 'String'), movieNo, 99999);
        plotType = get(main.dataDisplayType, 'Value');
        xmin = max(str2double(get(main.xAxisLimitLow, 'String')), 1);
        xmax = min(str2double(get(main.xAxisLimitHigh, 'String')), size(dataTraces,3));
        channelsToUse = getCluster512(centerChannel);

        
        if plotType == 2 %subtract mean
            subtractionVector = mean(dataTraces(:, channelsToUse, xmin:xmax), 1);
            dataToCluster = zeros(size(dataTraces, 1), length(channelsToUse), xmax-xmin+1);
            for i = 1:size(dataTraces,1)
                dataToCluster(i, :, :) = dataTraces(i, channelsToUse, xmin:xmax) - subtractionVector;
            end
        else
            dataToCluster = dataTraces(:, channelsToUse, xmin:xmax);
            if plotType == 3 %subtract estimated artifact
                warnh = warndlg('No artifact estimate has been calculated.  Using raw data.');
                uiwait(warnh)
            end
        end
        
        chosenElec = chooseElecForManualClusterGui(main, channelsToUse);
        
        if isnan(patternNo)
            main.tempAnalysisNoSave = manualCluster2D(reshape(dataToCluster(:,channelsToUse == chosenElec,:),...
                size(dataToCluster, 1), []));
        else
            main.tempAnalysis{patternNo, movieNo} = manualCluster2D(reshape(dataToCluster(:,channelsToUse == chosenElec,:),...
                size(dataToCluster, 1), []));
        end
       
        [main, elecResp] = refreshAll(main);
    end


    function resetForNewData(hObject, eventdata) %#ok<INUSD> %resets hand clustering in explore mode when switching to new dataset
        main.tempAnalysis = [];
        main.tempAnalysisNoSave = [];
    end

    function fitToErf(hObject, eventdata) %#ok<INUSD>
        filePath = get(main.filePath, 'String');
        fileName = get(main.fileName, 'String');      
        temp = load([filePath fileName]);
        elecResp = temp.elecResp;
        
        if isfield(elecResp, 'stimInfo') %if elecResp is loaded up
            nMovies = length(elecResp.stimInfo.movieNos);
             
            %checks to make sure values in successRates are correct
            for i = 1:nMovies
                mNum = elecResp.stimInfo.movieNos(i);
                if elecResp.stimInfo.nPulses(i) == 0 %value hasn't been set yet
                    dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo, mNum, 99999);
                    elecResp.stimInfo.nPulses(i) = size(dataTraces, 1);
                elseif ~isempty(elecResp.analysis.type{i}) && elecResp.stimInfo.nPulses(i) ~= length(elecResp.analysis.latencies{i})
                    warndlg('number of pulses listed in elecResp.stimInfo.nPulses(i) ~= length of elecResp.analysis.latencies{i}')
                    dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo, mNum, 5000);
                    elecResp.stimInfo.nPulses(i) = size(dataTraces, 1);
                    if elecResp.stimInfo.nPulses(i) ~= length(elecResp.analysis.latencies{i})
                        errordlg('number of pulses in preprocessed data does not match length of elecResp.analysis.latencies{i}')
                    end
                end
                if ~isempty(elecResp.analysis.type{i})
                    successRate = sum(elecResp.analysis.latencies{i}~=0)/elecResp.stimInfo.nPulses(i);
                    if elecResp.analysis.successRates(i) ~= successRate
                        disp(['warning: success rate for movie ' num2str(mNum) ' had to be corrected'])
                        elecResp.analysis.successRates(i) = successRate;
                    end
                end
            end
           
            
            data = zeros(2, nMovies);
            data(1,:) = elecResp.stimInfo.stimAmps;
            data(2,:) = elecResp.analysis.successRates;
            data(3,:) = elecResp.stimInfo.nPulses;
            lockedAmps = elecResp.analysis.finalized;
            for i = length(elecResp.stimInfo.stimAmps): -1: 1
                if isempty(elecResp.analysis.type{i})
                    data(:,i) = [];
                    lockedAmps(i) = [];
                end
            end
        else
            warnH = warndlg('need to load up an elecResp file first!');
            uiwait(warnH)
        end
        
        % linear-based
        data(1,:) = abs(data(1,:));
        [elecResp.analysis.erfParams, ~, elecResp.analysis.erfErr] = erfFitter(data, 2, -1, 'makePlot', 1, 'lockedAmps', lockedAmps);
        elecResp.analysis.threshold = -elecResp.analysis.erfParams(2)/elecResp.analysis.erfParams(1);
                
        % determines the minumum possible threshold for incomplete curves
%         if ~any(data(2,:) > 0.7)
%             minThreshData = data;
%             for ii = 1:5
%                 x = size(minThreshData,2);
%                 minThreshData(1,x+1) = minThreshData(1,x)*1.1;
%                 minThreshData(2,x+1) = 1;
%                 minThreshData(3,x+1) = mean(data(3,:));
%                 lockedAmps(x+1) = 0;
%             end
%             
%             
%             [~] = erfFitter(minThreshData, 2, -1, 'makePlot', 1, 'lockedAmps', lockedAmps);
%             hold on
%             text(min(data(1,:)), 0.9, 'MINIMUM POSSIBLE THRESH')
%             hold off          
%         end

        
        %standard deviations: don't redo if analysis hasn't changed
        if ~elecResp.analysis.erfCurrent
            elecResp.analysis.threshStd = bootstrapThresh(elecResp, main.bootstrapReps);
            elecResp.analysis.details.bootstrapReps = main.bootstrapReps;
        end
        
        elecResp.analysis.erfCurrent = 1;
        
        filePath = get(main.filePath, 'String');
        fileName = get(main.fileName, 'String');
        save([filePath filesep fileName], 'elecResp')
    end


    function altPlots(~, ~)
        altClusterPlots(main)
    end





    function switchMode(hObject, eventdata) %#ok<INUSD>
        
        if get(main.analyzeModeButton, 'Value')
            set(main.filePathStaticText, 'String', 'path to elecResp .mat file')
            set(main.fileNameStaticText, 'String', 'name of elecResp .mat file')
            set(main.residLimStaticText, 'String', 'analysis region limits')
            set(main.tempOffsetStaticText, 'String', 'template offset limits')
            set(main.tempStepSizeStaticText, 'String', 'step size:')
            set(main.excludeRegStaticText', 'Visible', 'on')
            
            set(main.fileName, 'visible', 'on', 'enable', 'on')
            set(main.eiPath, 'visible', 'off', 'enable', 'off')
            set(main.centerElec, 'visible', 'off', 'enable', 'off')
            set(main.residualLimitLow, 'enable', 'on', 'visible', 'on')
            set(main.residualLimitHigh, 'enable', 'on', 'Visible', 'on')
            set(main.offsetLimitLow, 'enable', 'on', 'Visible', 'on')
            set(main.offsetLimitHigh, 'enable', 'on', 'Visible', 'on')
            set(main.offsetStep, 'enable', 'on', 'Visible', 'on')
            set(main.spikeMinExclusionReg, 'enable', 'on', 'Visible', 'on')
            
            set(main.reanalyizeStaticText, 'visible', 'on')
            set(main.reanalyzeType, 'enable', 'on', 'visible', 'on')
            set(main.clearFlagsCheckBox, 'enable', 'on', 'visible', 'on')
            set(main.noRefitCheckBox, 'enable', 'on', 'visible', 'on')
            set(main.multipleMoviesCheckBox, 'enable', 'on', 'visible', 'on')
            set(main.reanalyzeButton, 'enable', 'on', 'visible', 'on')
            set(main.fitErfButton, 'enable', 'on')
            set(main.addActiveNeuron, 'String', 'add active neuron')
            set(main.removeFPButton, 'enable', 'on')
            set(main.removeFPButton, 'enable', 'on')
            set(main.addMissButton, 'enable', 'on')
            set(main.lockLatenciesButton, 'enable', 'on')
            set(main.clearAnalysisButton, 'enable', 'on')
            set(main.clearRemainingAnalysisButton, 'enable', 'on')
            set(main.clearPrecedingAnalysisButton, 'enable', 'on')
            
            set(main.exploreAnalyzeButton, 'visible', 'off', 'enable', 'off')
            
            set(main.patternNoStaticText, 'Visible', 'off')
            set(main.patternNo, 'Visible', 'off', 'enable', 'off')
            set(main.movieNoStaticText, 'Position', [10 110 130 15], 'String', 'movie number')
            set(main.movieNo, 'Position', [150 110 50 25])
            set(main.downButton, 'Position', [210 110 24 24])
            set(main.upButton, 'Position', [240 110 24 24])
            set(main.increaseOccButton, 'visible', 'on', 'enable', 'on')
            
            set(main.multiElecRespButton, 'visible', 'on', 'enable', 'on')
            
            set(main.fixedLimsCB, 'enable', 'on', 'visible', 'on')
            set(main.fixedLimsCB2, 'visible', 'on')
            switchLimMode()
            
            set(main.exploreOtherElecsCheckBox, 'visible', 'on', 'enable', 'on')
            set(main.bottomSpacer1, 'visible', 'on')
            set(main.bottomSpacer2, 'visible', 'on')
            if get(main.exploreOtherElecsCheckBox, 'Value')
                set(main.centerExploreNoSpikes, 'Visible', 'on', 'enable', 'on')
                set(main.subtractSpikesCheckBox, 'Visible', 'on', 'enable', 'on')
                set(main.excludeSpikeTracesCheckBox, 'Visible', 'on', 'enable', 'on');
                set(main.bottomPanelText, 'Visible', 'on')
            end
            
            
        else
            elecRespDummy.dummyField = [];
            elecResp = elecRespDummy;
            set(main.filePathStaticText, 'String', 'path to preprocessed data')
            set(main.fileNameStaticText, 'String', 'path to ei file')
            set(main.residLimStaticText, 'String', 'center electrode')
            set(main.tempOffsetStaticText, 'String', '')
            set(main.tempStepSizeStaticText, 'String', '')
            set(main.excludeRegStaticText', 'Visible', 'off')
            
            set(main.fileName, 'visible', 'off', 'enable', 'off')
            set(main.eiPath, 'visible', 'on', 'enable', 'on')
            set(main.centerElec, 'visible', 'on', 'enable', 'on')
            set(main.residualLimitLow, 'enable', 'off', 'visible', 'off')
            set(main.residualLimitHigh, 'enable', 'off', 'Visible', 'off')
            set(main.offsetLimitLow, 'enable', 'off', 'Visible', 'off')
            set(main.offsetLimitHigh, 'enable', 'off', 'Visible', 'off')
            set(main.offsetStep, 'enable', 'off', 'Visible', 'off')
            set(main.spikeMinExclusionReg, 'enable', 'off', 'Visible', 'off')
            
            set(main.reanalyizeStaticText, 'visible', 'off')
            set(main.reanalyzeType, 'enable', 'off', 'visible', 'off')
            set(main.clearFlagsCheckBox, 'enable', 'off', 'visible', 'off')
            set(main.noRefitCheckBox, 'enable', 'off', 'visible', 'off')
            set(main.multipleMoviesCheckBox, 'enable', 'off', 'visible', 'off')
            set(main.reanalyzeButton, 'enable', 'off', 'visible', 'off')
            set(main.fitErfButton, 'enable', 'off')
            set(main.addActiveNeuron, 'String', 'create new ei')
            set(main.removeFPButton, 'enable', 'off')
            set(main.removeFPButton, 'enable', 'off')
            set(main.addMissButton, 'enable', 'off')
            set(main.lockLatenciesButton, 'enable', 'off')
            set(main.clearAnalysisButton, 'enable', 'off')
            set(main.clearRemainingAnalysisButton, 'enable', 'off')
            set(main.clearPrecedingAnalysisButton, 'enable', 'off')
            
            set(main.exploreAnalyzeButton, 'visible', 'on', 'enable', 'on')
            
            set(main.patternNoStaticText, 'Visible', 'on')
            set(main.patternNo, 'Visible', 'on', 'enable', 'on')
            set(main.movieNoStaticText, 'Position', [205 110 35 15], 'String', 'movie')
            set(main.movieNo, 'Position', [245 110 50 25])
            set(main.downButton, 'Position', [300 110 24 24])
            set(main.upButton, 'Position', [325 110 24 24])
            set(main.increaseOccButton, 'visible', 'off', 'enable', 'off')
            
            set(main.centerExploreNoSpikes, 'Visible', 'off', 'enable', 'off')
            set(main.subtractSpikesCheckBox, 'Visible', 'off', 'enable', 'off')
            set(main.excludeSpikeTracesCheckBox, 'Visible', 'off', 'enable', 'off');
            
            set(main.bottomPanelText, 'Visible', 'off')
            set(main.bottomSpacer1, 'visible', 'off')
            set(main.bottomSpacer2, 'visible', 'off')
            set(main.exploreOtherElecsCheckBox, 'visible', 'off', 'enable', 'off')
            
            set(main.multiElecRespButton, 'visible', 'off', 'enable', 'off')
            
            set(main.fixedLimsCB, 'enable', 'off', 'visible', 'off')
            set(main.fixedLimsCB2, 'visible', 'off')
            switchLimMode()


        end
    end

    function switchLimMode(~,~)
        if get(main.fixedLimsCB, 'value')
            %set(main.xAxisLimitLow, 'enable', 'off')
            %set(main.xAxisLimitHigh, 'enable', 'off')
            set(main.residualLimitLow, 'enable', 'off')
            set(main.residualLimitHigh, 'enable', 'off')
            set(main.offsetLimitLow, 'enable', 'off')
            set(main.offsetLimitHigh, 'enable', 'off')
            set(main.offsetStep, 'enable', 'off')
            
            filePath = get(main.filePath, 'String');
            fileName = get(main.fileName, 'String');
            temp = load([filePath fileName]);
            
            movieIndex = find(temp.elecResp.stimInfo.movieNos == str2double(get(main.movieNo, 'String')));
            if ~isempty(temp.elecResp.analysis.details.residCalcWindow{movieIndex})
                
                spikeMinWindow = convertTempOffsetWindow(temp.elecResp, movieIndex);                
                
                %set(main.xAxisLimitLow, 'string', num2str(temp.elecResp.analysis.details.residCalcWindow{movieIndex}(1)-10))
                %set(main.xAxisLimitHigh, 'string', num2str(temp.elecResp.analysis.details.residCalcWindow{movieIndex}(2)+5))
                set(main.residualLimitLow, 'string', num2str(temp.elecResp.analysis.details.residCalcWindow{movieIndex}(1)))
                set(main.residualLimitHigh, 'string', num2str(temp.elecResp.analysis.details.residCalcWindow{movieIndex}(2)))
                set(main.offsetLimitLow, 'string', num2str(spikeMinWindow(1)))
                set(main.offsetLimitHigh, 'string', num2str(spikeMinWindow(2)))
                set(main.offsetStep, 'string', num2str(temp.elecResp.analysis.details.tempOffsetStep{movieIndex}))
            end

        else
            %set(main.xAxisLimitLow, 'enable', 'on')
            %set(main.xAxisLimitHigh, 'enable', 'on')
            set(main.residualLimitLow, 'enable', 'on')
            set(main.residualLimitHigh, 'enable', 'on')
            set(main.offsetLimitLow, 'enable', 'on')
            set(main.offsetLimitHigh, 'enable', 'on')
            set(main.offsetStep, 'enable', 'on')
        end
    end

    function otherElecsSwitch(hObject, eventdata) %#ok<INUSD>
        if get(main.exploreOtherElecsCheckBox, 'Value')
            set(main.bottomPanelText, 'visible', 'on')
            set(main.centerExploreNoSpikes, 'visible', 'on', 'enable', 'on')
            set(main.subtractSpikesCheckBox, 'visible', 'on', 'enable', 'on')
            set(main.excludeSpikeTracesCheckBox, 'Visible', 'on', 'enable', 'on');
        else
            set(main.bottomPanelText, 'visible', 'off')
            set(main.centerExploreNoSpikes, 'visible', 'off', 'enable', 'off')
            set(main.subtractSpikesCheckBox, 'visible', 'off', 'enable', 'off')
            set(main.excludeSpikeTracesCheckBox, 'Visible', 'off', 'enable', 'off');
        end
    end
end