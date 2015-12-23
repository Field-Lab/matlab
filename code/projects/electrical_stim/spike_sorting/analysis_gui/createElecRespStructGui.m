function createElecRespStructGui(main)
% makes a gui form that can be used to generate an elecResp.mat filemin

%% parameters

interval = 50; %length of time before stimulus is updated, in microseconds
elecRespH.editTextFont = 'Courier';
elecRespH.displayTextFontSize = 9;

%% preparations

filePath = get(main.filePath, 'String');
stringIndeces = strfind(filePath, '/');
eiPath = get(main.eiPath, 'String');
stringIndecesEi = strfind(eiPath, '/');
elecRespH.autoMovieVal = 1;

%%

elecRespH.gui = figure('position', [300 400 600 565], 'Visible', 'off');

uicontrol(elecRespH.gui, 'Style', 'text', 'Position', [10 545 580 15],...
    'String', ' Required', 'HorizontalAlignment', 'left', 'Fontsize', 14,...
    'BackgroundColor', [0.8 0 0], 'ForegroundColor', 'white');
elecRespH.requiredPanel = uipanel(elecRespH.gui, 'Units', 'pixels', 'Position', [10 280 580 275],...
    'BackgroundColor', 'white', 'Fontsize', 10, 'BorderType', 'line', 'HighlightColor', 'black',...
    'BorderWidth', 0);

uicontrol(elecRespH.gui, 'Style', 'text', 'Position', [10 250 580 15],...
    'String', ' Optional', 'HorizontalAlignment', 'left', 'Fontsize', 14,...
    'BackgroundColor', [0.8 0 0], 'ForegroundColor', 'white');
elecRespH.optionalPanel = uipanel(elecRespH.gui, 'Units', 'pixels', 'Position', [10 60 580 190],...
    'BackgroundColor', 'white', 'Fontsize', 10, 'BorderType', 'line', 'HighlightColor', 'black',...
    'BorderWidth', 0);

% pushbuttons
elecRespH.okButton = uicontrol(elecRespH.gui,  'Style', 'pushbutton', 'String', 'create elecResp file',...
        'Position', [100 20 150 20], 'Callback', @saveElecResp, 'FontWeight', 'bold');

elecRespH.cancelButton = uicontrol(elecRespH.gui,  'Style', 'pushbutton', 'String', 'cancel',...
        'Position', [350 20 150 20], 'Callback', @cancelFun, 'FontWeight', 'bold');


%% static text (required panel)

% paths
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [10 235 200 15],...
    'String', 'path to preprocessed data (directory): ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [10 210 200 15],...
    'String', 'path to VISION analysis files: ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [10 185 200 15],...
    'String', 'base name of VISION analysis files: ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [10 172 200 15],...
    'String', '(e.g. data000) ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', 7);
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [10 160 200 15],...
    'String', 'destination for saved files: ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);

% pattern, movie numbers
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [10 130 200 15],...
    'String', 'pattern number: ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [335 130 105 15],...
    'String', 'movie numbers:', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [360 80 30 15],...
    'String', 'first ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [425 80 30 15],...
    'String', 'last ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [495 80 40 15],...
    'String', 'interval ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);

% neurons
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [10 40 200 15],...
    'String', 'target neuron: ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [255 40 135 15],...
    'String', 'other active neurons: ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);

% electrodes
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [10 10 200 15],...
    'String', 'main recording electrode: ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.requiredPanel, 'Style', 'text', 'Position', [255 10 135 15],...
    'String', 'other electrodes to use: ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);


%% edit text boxes (required panel)

% paths
elecRespH.dataPathText = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', filePath,...
    'Position', [215 235 355 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.analysisPathText = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', eiPath(1:stringIndecesEi(end)),...
    'Position', [215 210 355 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.analysisBaseNameText = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', eiPath(stringIndecesEi(end)+1:strfind(eiPath, '.')-1),...
    'Position', [215 185 355 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.savePathText = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', filePath,...
    'Position', [215 160 355 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);

% pattern, movie numbers
elecRespH.patternNoText = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', get(main.patternNo, 'String'),...
    'Position', [215 130 130 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont,  'Callback', @enterNeuronOrPattern);
elecRespH.movieFirstText = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '1',...
    'Position', [395 80 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont, 'enable', 'off');
elecRespH.movieLastText = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '1',...
    'Position', [460 80 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont, 'enable', 'off');
elecRespH.movieIntText = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '1',...
    'Position', [540 80 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont, 'enable', 'off');

elecRespH.autoMovieCheck = uicontrol(elecRespH.requiredPanel, 'Style', 'checkbox', 'String', 'auto',...
    'Position', [450 130 100 25], 'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize,...
    'Value', 1, 'Callback', @autoMovieCB);

elecRespH.autoButtonGroup = uibuttongroup('Parent', elecRespH.requiredPanel, 'Units', 'pixels',...
    'Position', [350 110 230 20], 'BorderType', 'none', 'BackgroundColor', 'white');
elecRespH.allButton = uicontrol(elecRespH.autoButtonGroup, 'Style', 'radiobutton', 'String', 'all', 'Position', [0 0 40 20]);
elecRespH.oddButton = uicontrol(elecRespH.autoButtonGroup, 'Style', 'radiobutton', 'String', 'odds', 'Position', [40 0 60 20]);
elecRespH.evenButton = uicontrol(elecRespH.autoButtonGroup, 'Style', 'radiobutton', 'String', 'evens', 'Position', [105 0 60 20]);
elecRespH.bothButton = uicontrol(elecRespH.autoButtonGroup, 'Style', 'radiobutton', 'String', 'both', 'Position', [170 0 60 20]);

% neurons
elecRespH.mainNeuron = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '1',...
    'Position', [215 40 40 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont, 'Callback', @enterNeuronOrPattern);
elecRespH.activeNeurons{1} = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '',...
   'Position', [395 40 40 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.activeNeurons{2} = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '',...
   'Position', [440 40 40 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.activeNeurons{3} = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '',...
   'Position', [485 40 40 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.activeNeurons{4} = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '',...
   'Position', [530 40 40 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);

elecRespH.externalEiCheck = uicontrol(elecRespH.requiredPanel, 'Style', 'checkbox', 'String', 'external ei file',...
    'Position', [215 70 100 25], 'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize,...
    'Value', 0);

% electrodes
elecRespH.mainElec = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', get(main.centerElec, 'String'),...
    'Position', [215 10 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.otherElecs{1} = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '',...
    'Position', [395 10 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.otherElecs{2} = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '',...
    'Position', [430 10 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.otherElecs{3} = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '',...
    'Position', [465 10 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.otherElecs{4} = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '',...
    'Position', [500 10 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.otherElecs{5} = uicontrol(elecRespH.requiredPanel, 'Style', 'edit', 'String', '',...
    'Position', [535 10 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);

%% static text (optional)

uicontrol(elecRespH.optionalPanel, 'Style', 'text', 'Position', [10 160 200 15],...
    'String', 'primary electrode: ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.optionalPanel, 'Style', 'text', 'Position', [10 130 200 15],...
    'String', 'experiment name: ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.optionalPanel, 'Style', 'text', 'Position', [10 100 200 15],...
    'String', 'dataset short name: ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);

uicontrol(elecRespH.optionalPanel, 'Style', 'text', 'Position', [10 70 200 15],...
    'String', 'path to ttx data (directory): ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);

% artifact movie numbers
uicontrol(elecRespH.optionalPanel, 'Style', 'text', 'Position', [10 40 200 15],...
    'String', 'ttx movie numbers:', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.optionalPanel, 'Style', 'text', 'Position', [230 40 30 15],...
    'String', 'first ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.optionalPanel, 'Style', 'text', 'Position', [300 40 30 15],...
    'String', 'last ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);
uicontrol(elecRespH.optionalPanel, 'Style', 'text', 'Position', [380 40 40 15],...
    'String', 'interval ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white', 'Fontsize', elecRespH.displayTextFontSize);

uicontrol(elecRespH.optionalPanel, 'Style', 'text', 'Position', [10 10 200 15],...
    'String', 'sample rate (Hz): ', 'HorizontalAlignment', 'right',...
    'BackgroundColor', 'white');


%% edit text (optional)
elecRespH.pElecText = uicontrol(elecRespH.optionalPanel, 'Style', 'edit', 'String', '',...
    'Position', [215 160 355 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.experimentNameText = uicontrol(elecRespH.optionalPanel, 'Style', 'edit', 'String', filePath(stringIndeces(end-2)+1:stringIndeces(end-1)-1),...
    'Position', [215 130 355 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.shortNameText = uicontrol(elecRespH.optionalPanel, 'Style', 'edit', 'String', [filePath(stringIndeces(end-2)+1:stringIndeces(end)-1) '_nX_p' get(main.patternNo, 'String')],...
    'Position', [215 100 355 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);


elecRespH.artifactPathText = uicontrol(elecRespH.optionalPanel, 'Style', 'edit', 'String', '',...
    'Position', [215 70 355 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);

% artifact movie numbers
elecRespH.artMovieFirstText = uicontrol(elecRespH.optionalPanel, 'Style', 'edit', 'String', '',...
    'Position', [265 40 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.artMovieLastText = uicontrol(elecRespH.optionalPanel, 'Style', 'edit', 'String', '',...
    'Position', [335 40 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);
elecRespH.artMovieIntText = uicontrol(elecRespH.optionalPanel, 'Style', 'edit', 'String', '',...
    'Position', [425 40 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);

% sample rate
elecRespH.sampleRateText = uicontrol(elecRespH.optionalPanel, 'Style', 'edit', 'String', '20000',...
    'Position', [215 10 60 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white', 'FontName', elecRespH.editTextFont);



%% initialize gui

set(elecRespH.gui, 'Visible', 'on')

uiwait(elecRespH.gui)

%% callbacks

    function enterNeuronOrPattern(hObject, eventdata) %#ok<INUSD>'
        set(elecRespH.shortNameText, 'String',...
            [filePath(stringIndeces(end-2)+1:stringIndeces(end)-1)...
            '_n' get(elecRespH.mainNeuron, 'String') '_p' get(elecRespH.patternNoText, 'String')])    
    end

    function autoMovieCB(hObject, eventdata) %#ok<INUSD>
        elecRespH.autoMovieVal = get(elecRespH.autoMovieCheck, 'Value');
        if elecRespH.autoMovieVal
            set(elecRespH.movieFirstText, 'enable', 'off')
            set(elecRespH.movieLastText, 'enable', 'off')
            set(elecRespH.movieIntText, 'enable', 'off')
            
            set(elecRespH.allButton, 'enable', 'on')
            set(elecRespH.oddButton, 'enable', 'on')
            set(elecRespH.evenButton, 'enable', 'on')
            set(elecRespH.bothButton, 'enable', 'on')
        else
            set(elecRespH.movieFirstText, 'enable', 'on')
            set(elecRespH.movieLastText, 'enable', 'on')
            set(elecRespH.movieIntText, 'enable', 'on')
            
            set(elecRespH.allButton, 'enable', 'off')
            set(elecRespH.oddButton, 'enable', 'off')
            set(elecRespH.evenButton, 'enable', 'off')
            set(elecRespH.bothButton, 'enable', 'off')
        end
        
    end

    function saveElecResp(hObject, eventdata) %#ok<INUSD>
        
        elecRespInfo.dataPath = get(elecRespH.dataPathText, 'String');
        

        elecRespInfo.analysisPath =     get(elecRespH.analysisPathText, 'String');
        elecRespInfo.analysisBaseName = get(elecRespH.analysisBaseNameText, 'String');
        elecRespInfo.savePath =         get(elecRespH.savePathText, 'String');

        elecRespInfo.patternNo =  str2double(get(elecRespH.patternNoText, 'String'));
        if isnan(elecRespInfo.patternNo) %in case pattern number includes other characters
            elecRespInfo.patternNo =  get(elecRespH.patternNoText, 'String');
        end
        elecRespInfo.movieFirst = str2double(get(elecRespH.movieFirstText, 'String'));
        elecRespInfo.movieLast =  str2double(get(elecRespH.movieLastText, 'String'));
        elecRespInfo.movieInt =   str2double(get(elecRespH.movieIntText, 'String'));

        elecRespInfo.mainNeuron =       str2double(get(elecRespH.mainNeuron, 'String'));
        elecRespInfo.activeNeurons{1} = str2double(get(elecRespH.activeNeurons{1}, 'String'));
        elecRespInfo.activeNeurons{2} = str2double(get(elecRespH.activeNeurons{2}, 'String'));
        elecRespInfo.activeNeurons{3} = str2double(get(elecRespH.activeNeurons{3}, 'String'));
        elecRespInfo.activeNeurons{4} = str2double(get(elecRespH.activeNeurons{4}, 'String'));
        elecRespInfo.externalEi = get(elecRespH.externalEiCheck, 'Value');

        elecRespInfo.mainElec =      str2double(get(elecRespH.mainElec, 'String'));
        elecRespInfo.otherElecs{1} = str2double(get(elecRespH.otherElecs{1}, 'String'));
        elecRespInfo.otherElecs{2} = str2double(get(elecRespH.otherElecs{2}, 'String'));
        elecRespInfo.otherElecs{3} = str2double(get(elecRespH.otherElecs{3}, 'String'));
        elecRespInfo.otherElecs{4} = str2double(get(elecRespH.otherElecs{4}, 'String'));
        elecRespInfo.otherElecs{5} = str2double(get(elecRespH.otherElecs{5}, 'String'));
        
        elecRespInfo.pElec =          str2double(get(elecRespH.pElecText, 'String'));
        elecRespInfo.experimentName = get(elecRespH.experimentNameText, 'String');
        elecRespInfo.shortName =      get(elecRespH.shortNameText, 'String');
        elecRespInfo.artifactPath =   get(elecRespH.artifactPathText, 'String');
        elecRespInfo.artMovieFirst =  str2double(get(elecRespH.artMovieFirstText, 'String'));
        elecRespInfo.artMovieLast =   str2double(get(elecRespH.artMovieLastText, 'String'));
        elecRespInfo.artMovieInt =    str2double(get(elecRespH.artMovieIntText, 'String'));
        elecRespInfo.sampleRate =     str2double(get(elecRespH.sampleRateText, 'String'));
        
        elecRespInfo.autoMovie = elecRespH.autoMovieVal;
        
        if elecRespInfo.autoMovie
            bothAutoTypes = false;
            selectedAuto = get(elecRespH.autoButtonGroup, 'SelectedObject');
            if selectedAuto == elecRespH.allButton
                elecRespInfo.autoType = 'all';
            elseif selectedAuto == elecRespH.oddButton
                elecRespInfo.autoType = 'odds';
            elseif selectedAuto == elecRespH.evenButton%selectedAuto == elecRespH.evenButton
                elecRespInfo.autoType = 'evens';
            else
                bothAutoTypes = true;
            end
        end
 
        if isnumeric(elecRespInfo.patternNo)
            patternNoString = num2str(elecRespInfo.patternNo);
        else
            patternNoString = elecRespInfo.patternNo;
        end
        
        saveName = ['elecResp_n' num2str(elecRespInfo.mainNeuron) '_p' patternNoString '.mat'];
        saveNameOdds = ['elecResp_n' num2str(elecRespInfo.mainNeuron) '_p' patternNoString '_w50.mat'];
        saveNameEvens = ['elecResp_n' num2str(elecRespInfo.mainNeuron) '_p' patternNoString '_w100.mat'];
        
        
        if elecRespInfo.autoMovie && bothAutoTypes
            elecRespInfo.autoType = 'odds';
            elecResp = createElecRespStruct(elecRespInfo);
            uisave('elecResp', saveNameOdds)
            
            elecRespInfo.autoType = 'evens';
            elecResp = createElecRespStruct(elecRespInfo);
            uisave('elecResp', saveNameEvens)
        else
            elecResp = createElecRespStruct(elecRespInfo);
            
            if elecRespInfo.autoMovie && strcmp(elecRespInfo.autoType, 'odds')
                uisave('elecResp', saveNameOdds)
            elseif elecRespInfo.autoMovie && strcmp(elecRespInfo.autoType, 'evens')
                uisave('elecResp', saveNameEvens)
            else
                uisave('elecResp', saveName)
            end
        end
        
%         if isempty(elecResp)
%             return
%         end
        
        if ~exist(elecResp.names.savePath, 'file')
            mkdir(elecResp.names.savePath)
        end

        close(elecRespH.gui)
    end


    function cancelFun(hObject, eventdata) %#ok<INUSD,INUSD>
        close(elecRespH.gui)
    end

end