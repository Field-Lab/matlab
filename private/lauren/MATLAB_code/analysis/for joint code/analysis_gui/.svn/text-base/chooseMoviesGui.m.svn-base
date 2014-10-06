function selectedMovies = chooseMoviesGui(elecResp)

selectedMovies = [];
movieNos = elecResp.stimInfo.movieNos;
nMovies = length(movieNos);

listLength = ceil(nMovies/2);

chooseMoviesGui = figure('position', [500 500 240 170+20*listLength], 'color', 'white', 'Toolbar', 'none',...
    'Menubar', 'none', 'visible', 'off');

uicontrol(chooseMoviesGui, 'Style', 'text', 'String', 'Choose movie numbers to analyze',...
    'Position', [20, 130+20*listLength, 200, 20], 'BackgroundColor', 'white', 'HorizontalAlignment', 'center')

movieCheckBox = cell(1, nMovies);
for i = 1:listLength
    movieCheckBox{i} = uicontrol(chooseMoviesGui, 'Style', 'checkbox', 'Units', 'pixels',...
        'Position', [50 100+20*listLength-(i-1)*20 90 20], 'String', num2str(movieNos(i)), 'Value', 0);
    if i ~= listLength || ~mod(nMovies, 2)
        movieCheckBox{listLength+i} = uicontrol(chooseMoviesGui, 'Style', 'checkbox', 'Units', 'pixels',...
            'Position', [150 100+20*listLength-(i-1)*20 90 20], 'String', num2str(movieNos(listLength+i)), 'Value', 0);
    end
end

uicontrol(chooseMoviesGui, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [20 85 200 20],...
    'String', 'check all', 'Callback', @checkAll)
uicontrol(chooseMoviesGui, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [20 60 200 20],...
    'String', 'uncheck all', 'Callback', @uncheckAll)
uicontrol(chooseMoviesGui, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [20 20 200 20],...
    'String', 'accept', 'Callback', @acceptMovies)

%% initialize

set(chooseMoviesGui, 'visible', 'on')

uiwait(chooseMoviesGui)

%% callbacks

    function checkAll(hObject, eventdata) %#ok<INUSD>
        for ii = 1:nMovies
            set(movieCheckBox{ii}, 'Value', 1)
        end
    end

    function uncheckAll(hObject, eventdata) %#ok<INUSD>
        for ii = 1:nMovies
            set(movieCheckBox{ii}, 'Value', 0)
        end
    end

    function acceptMovies(hObject, eventdata) %#ok<INUSD>
        selectedMovies = [];
        for ii = 1:nMovies
            if get(movieCheckBox{ii}, 'Value')
                selectedMovies = [selectedMovies movieNos(ii)];
            end
        end
        close(chooseMoviesGui)
    end


end