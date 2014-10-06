function [direction pts eCoords] = traceAxonDir(ei)

[xCoords yCoords] = getElectrodeCoords61();
nElecs = length(xCoords);

%scale x coordinates to reflect actual physical array
scaleFactor = 2/xCoords(1); %distance of electrode 1 from center should be 2 
%(where 2 = physical between vertical or horizontal pair of electrodes)
xCoords = xCoords*scaleFactor;

%package for function return
eCoords.x = xCoords;
eCoords.y = yCoords;

chooseAxonDir.gui = figure('position', [100 150 700 700], 'Toolbar', 'none', 'Menubar', 'none', 'Color', 'white',...
    'visible', 'off');

uicontrol(chooseAxonDir.gui,  'Style', 'pushbutton',...
    'String', 'finalize selection', 'Position', [400 20 250 20], 'Callback', @finalize);

uicontrol(chooseAxonDir.gui,  'Style', 'pushbutton',...
    'String', 'choose endpoints', 'Position', [50 20 250 20], 'Callback', @newClicks);

aEI = axes('parent', chooseAxonDir.gui, 'units', 'pixels', 'position', [50 60 600 600]); axis equal; hold on
title(aEI, 'click soma side first')

refreshEIPlot()


chooseAxonDir.x_pt = [0 0];
chooseAxonDir.y_pt = [0 0];

uiwait(chooseAxonDir.gui)

    function refreshEIPlot()
        cla(aEI)
        %plot EI: amplitude is based on peak negative voltage
        maxEISig = max(max(-ei));
        mSc = 50; %scale factor for marker size
        for ii = 1:nElecs;
            if ~any([9 25 57]==ii)
            plot(aEI, xCoords(ii), yCoords(ii), 'o', 'markerSize', mSc*max(-ei(ii,:))/maxEISig,...
                'markerFaceColor', [0 0 1], 'markerEdgeColor', [0 0 1])
            
            %scale y values so that axon signals aren't tiny, but pure
            %noise isn't too big
            %plot(aEI, xCoords(ii)+linspace(-1, 1, size(ei,2)), yCoords(ii)+ei(ii,:)/max(abs(ei(ii,:))), '-', 'color', [0.5 0.5 0.5])
            plot(aEI, xCoords(ii)+linspace(-1, 1, size(ei,2)), yCoords(ii)+ei(ii,:)/(max(abs(ei(ii,:)))+0.05*maxEISig), '-', 'color', [0.5 0.5 0.5])
            end
        end
        set(aEI, 'xLim', [min(xCoords)-2 max(xCoords)+2], 'ylim', [min(yCoords)-2 max(yCoords)+2],...
            'dataAspectRatio', [1 1 1], 'xTick', [], 'yTick', [], 'box', 'on')
    end

    function newClicks(~,~)
        refreshEIPlot()
        %select points
        x_pt=[]; y_pt=[];
        while length(x_pt) < 2
            [a,b,key]=ginput(1);
            x_pt=[x_pt;a];y_pt=[y_pt;b]; %#ok<AGROW,AGROW>
            plot(x_pt, y_pt, 'ro');
        end
        plot(x_pt, y_pt, 'r-');
        chooseAxonDir.x_pt = x_pt;
        chooseAxonDir.y_pt = y_pt;
    end

    function finalize(~, ~)
        x_pt = chooseAxonDir.x_pt;
        y_pt = chooseAxonDir.y_pt;
        vect = [x_pt(2)-x_pt(1) y_pt(2)-y_pt(1)];
        direction = (180/pi)*atan2(vect(2), vect(1));
        
        pts.x = x_pt;
        pts.y = y_pt;
        
        plot(aEI, 0, 0, 'ro', 'markerfacecolor', [1 0 0])
        plot(aEI, [0 5*cosd(direction)], [0 5*sind(direction)], 'r-', 'linewidth', 2); drawnow
        pause(1)
        close(chooseAxonDir.gui)
    end
end

