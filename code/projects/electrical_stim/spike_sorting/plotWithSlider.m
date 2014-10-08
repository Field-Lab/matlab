function plotWithSlider(axesToPlot)
% axesToPlot(axesToPlot) plots a set of axes in a single figure, with a slider to scroll through
% which axes is visible at any given time
%
% arguments
%    axesToPlot: cell array of handles to axes
%
% 

H.sliderGui = figure('Visible', 'off', 'Position', [50 100 1200 600], 'Name',...
    'move slider to see other plots', 'Resize', 'off');

nAxes = length(axesToPlot);
H.axes = axesToPlot;

for i = 1:nAxes
    set(H.axes{i},'Parent', H.sliderGui, 'Visible', 'off')
    set(findall(H.axes{i}), 'Visible', 'off')
end

H.slider = uicontrol(H.sliderGui, 'Style', 'slider', 'Max', nAxes,...
    'Min', 1, 'Value', 1, 'SliderStep', [1/nAxes, 5/nAxes],...
    'units', 'normal','Position', [0 0 1 0.05],...
    'Callback',{@sliderCB, H.axes});


%% initialization

set(findall(H.axes{1}), 'Visible', 'on')
set(H.sliderGui, 'Visible', 'on')

%% callbacks

    function sliderCB(hObject, eventdata, axesHandles)
        for ii = 1:length(H.axes)
            sliderValue = round(get(hObject, 'Value'));
            if sliderValue == i
                obs = findall(H.axes{i});
                keyboard
                set(H.axes{i}, 'Visible', 'on')
                set(obs,'Visible','on');
            else
                obs = findall(H.axes{i});
                set(H.axes{i}, 'Visible', 'off')
                set(obs,'Visible','off');                
            end
        end
    end
uiwait

end