function checkMultipleOccurrences(main)


%% construct the gui

checkMultOccGui.figH = figure('position', [500 300 400 480], 'color', 'white', 'Toolbar', 'none',...
    'Menubar', 'none', 'Name', 'check for errors in multiple occurrences of pattern');

%static text
uicontrol(checkMultOccGui.figH, 'Style', 'text', 'String', 'elecResp base name',...
    'Position', [20 440 150 15], 'HorizontalAlignment', 'right', 'BackgroundColor', 'white');
uicontrol(checkMultOccGui.figH, 'Style', 'text', 'String', '(use * to indicate where variable part of name goes)',...
    'Position', [20 420 360 15], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white');

%edit text
checkMultOccGui.baseName = uicontrol(checkMultOccGui.figH, 'Style', 'edit', 'String', get(main.fileName, 'string'),...
    'Position', [180 440 200 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');


uicontrol(checkMultOccGui.figH, 'Style', 'text', 'String', 'variable part of elecResp name:',...
    'Position', [20 385 360 15], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
checkMultOccGui.variType = uicontrol(checkMultOccGui.figH, 'Style', 'popupmenu',...
    'String', {'series of integers', 'enter each manually'},...
    'Value', 1, 'Position', [100 350 200 30], 'BackgroundColor', [1 1 1], 'FontSize', 9, 'callback', @switchVariType);


%variable portion = series of integers
checkMultOccGui.intText1 =  uicontrol(checkMultOccGui.figH, 'Style', 'text', 'String', 'first integer',...
    'Position', [120 300 70 15], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
checkMultOccGui.intText2 = uicontrol(checkMultOccGui.figH, 'Style', 'text', 'String', 'last integer',...
    'Position', [120 270 70 15], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
checkMultOccGui.intText3 = uicontrol(checkMultOccGui.figH, 'Style', 'text', 'String', 'increment',...
    'Position', [120 220 70 15], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');

checkMultOccGui.startInt = uicontrol(checkMultOccGui.figH, 'Style', 'edit', 'String', '1',...
    'Position', [190 300 100 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'callback', @enteredStartInt);
checkMultOccGui.endInt = uicontrol(checkMultOccGui.figH, 'Style', 'edit', 'String', '10',...
    'Position', [190 270 100 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
checkMultOccGui.incInt = uicontrol(checkMultOccGui.figH, 'Style', 'edit', 'String', '1',...
    'Position', [190 220 100 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');



%manually entered variable portions

checkMultOccGui.variText = cell(1,20);
for ii = 1:10
    checkMultOccGui.variText{ii} = uicontrol(checkMultOccGui.figH, 'Style', 'edit', 'String', '',...
    'Position', [90 340-ii*25 100 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white',...
    'enable', 'off', 'visible', 'off', 'callback', @enteredVariString);
end
for ii = 11:20
    checkMultOccGui.variText{ii} = uicontrol(checkMultOccGui.figH, 'Style', 'edit', 'String', '',...
    'Position', [210 340+250-ii*25 100 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white',...
    'enable', 'off', 'visible', 'off', 'callback', @enteredVariString);
end

uicontrol(checkMultOccGui.figH, 'Style', 'text', 'String', 'example:',...
    'Position', [80 65 80 15], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
checkMultOccGui.exampleFile = uicontrol(checkMultOccGui.figH, 'Style', 'text', 'String', '',...
    'Position', [170 65 200 15], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');


uicontrol(checkMultOccGui.figH, 'Style', 'pushbutton', 'String', 'go!',...
    'Position', [20 20 170 25], 'Callback', @checkOccurrencesDummy)
uicontrol(checkMultOccGui.figH, 'Style', 'pushbutton', 'String', 'close',... 
    'Position', [210 20 170 25], 'Callback', @cancel)


%%

enteredStartInt()
uiwait(checkMultOccGui.figH)

    function checkOccurrencesDummy(hObject, eventdata) %#ok<INUSD>
        baseString = get(checkMultOccGui.baseName, 'string');
        elecRespList = cell(0);
        variType = get(checkMultOccGui.variType, 'value');
        
        if variType == 1 %integer sequence
            startInt = str2double(get(checkMultOccGui.startInt, 'string'));
            endInt = str2double(get(checkMultOccGui.endInt, 'string'));
            incInt = str2double(get(checkMultOccGui.incInt, 'string'));
            if isnan(startInt) || isnan(endInt) || isnan(incInt)
                warndlg('one of the values is not a valid integer')
                return
            end
            integers = startInt:incInt:endInt;
            for jj = 1:length(integers)
                tmpString = num2str(integers(jj));
                elecRespList{jj} = combineStrings(baseString, tmpString);
            end
        else %manual list
            kk = 0;
            for jj = 1:20
                tmpString = get(checkMultOccGui.variText{jj}, 'string');
                if ~isempty(tmpString)
                    kk = kk+1;
                    elecRespList{kk} = combineStrings(baseString, tmpString);
                end
            end
        end
        
        %make sure all listed elecResp files exist
        elecRespPath = get(main.filePath, 'String');
        for jj = 1:length(elecRespList)
            if ~exist([elecRespPath elecRespList{jj}], 'file')
                warndlg(['unable to find ' elecRespPath elecRespList{jj} ' and possibly others'])
                return
            end
        end
        
        checkOccurrencesPlotter(main, elecRespList)
    end

    function switchVariType(hObject, eventdata) %#ok<INUSD>
        variType = get(checkMultOccGui.variType, 'value');
        if variType == 1 %integer sequence
            set(checkMultOccGui.intText1, 'visible', 'on')
            set(checkMultOccGui.intText2, 'visible', 'on')
            set(checkMultOccGui.intText3, 'visible', 'on')
            set(checkMultOccGui.startInt, 'visible', 'on', 'enable', 'on')
            set(checkMultOccGui.endInt, 'visible', 'on', 'enable', 'on')
            set(checkMultOccGui.incInt, 'visible', 'on', 'enable', 'on')
            for jj = 1:20
                set(checkMultOccGui.variText{jj}, 'visible', 'off', 'enable', 'off')
            end
            enteredStartInt()

        else %manual list
            set(checkMultOccGui.intText1, 'visible', 'off')
            set(checkMultOccGui.intText2, 'visible', 'off')
            set(checkMultOccGui.intText3, 'visible', 'off')
            set(checkMultOccGui.startInt, 'visible', 'off', 'enable', 'off')
            set(checkMultOccGui.endInt, 'visible', 'off', 'enable', 'off')
            set(checkMultOccGui.incInt, 'visible', 'off', 'enable', 'off')
            for jj = 1:20
                set(checkMultOccGui.variText{jj}, 'visible', 'on', 'enable', 'on')
            end
            enteredVariString()
        end
    end

    function cancel(hObject, eventdata) %#ok<INUSD>
        close(checkMultOccGui.figH)
    end

    function enteredVariString(hObject, eventdata) %#ok<INUSD>
        tmp = get(checkMultOccGui.variText{1}, 'string');
        if ~isempty(tmp)
            baseString = get(checkMultOccGui.baseName, 'string');
            fullString = combineStrings(baseString, tmp);
            set(checkMultOccGui.exampleFile, 'string', fullString)
        else
            set(checkMultOccGui.exampleFile, 'string', '')
        end
    end

    function enteredStartInt(hObject, eventdata) %#ok<INUSD>
        tmp = get(checkMultOccGui.startInt, 'string');
        tmp = str2double(tmp);
        if ~isnan(tmp)
            baseString = get(checkMultOccGui.baseName, 'string');
            fullString = combineStrings(baseString, num2str(tmp));
            set(checkMultOccGui.exampleFile, 'string', fullString)
        else
            set(checkMultOccGui.exampleFile, 'string', '')
        end
    end



    function fullString = combineStrings(baseString, varString)
        astInd = strfind(baseString, '*');
        baseStringPre = baseString(1:astInd-1);
        baseStringPost = baseString(astInd+1:end);
        fullString = [baseStringPre varString baseStringPost];
    end
end