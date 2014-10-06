function [PCx PCy] = PCChooser(PCAScore)

% allows the user to choose which PC plot to use (for handclustering)
%
% arguments:
%     PCAScore must be a nxn array with n>=4 (only the first 4 PCs are used)
%
% returns 
%     PCx, PCy: the PCs identities of the PCs used in the chosen plot. PCx = PCy = 0 signifies that
%     the button "cluster all remaining" was pressed.
%
% written by Lauren Hruby 2008-09-08
% last edited 2008-09-09 


%% Construct the gui components

h.gui = figure('Visible', 'off', 'Position',[100 100 800 600], 'Name', 'Choose which PC plot gives the best cluster separation.', 'Resize', 'off');

h.p1h = axes('Parent', h.gui, 'Units', 'pixels', 'Position', [50 300 200 200]); %from left, from bottom, width, height
h.p2h = axes('Parent', h.gui, 'Units', 'pixels', 'Position', [300 300 200 200]);
h.p3h = axes('Parent', h.gui, 'Units', 'pixels', 'Position', [550 300 200 200]);
h.p4h = axes('Parent', h.gui, 'Units', 'pixels', 'Position', [50 50 200 200]);
h.p5h = axes('Parent', h.gui, 'Units', 'pixels', 'Position', [300 50 200 200]);
h.p6h = axes('Parent', h.gui, 'Units', 'pixels', 'Position', [550 50 200 200]);


uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'PC1 vs PC2', 'Position', [50 500 200 30], 'Callback', @PC12);
uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'PC1 vs PC3', 'Position', [300 500 200 30], 'Callback', @PC13);
uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'PC1 vs PC4', 'Position', [550 500 200 30], 'Callback', @PC14);
uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'PC2 vs PC3', 'Position', [50 250 200 30], 'Callback', @PC23);
uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'PC2 vs PC4', 'Position', [300 250 200 30], 'Callback', @PC24);
uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'PC3 vs PC4', 'Position', [550 250 200 30], 'Callback', @PC34);
uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'cluster all remaining', 'Position', [300 550 200 30], 'Callback', @clusterAll);

%% Initialization tasks

axes(h.p1h)
plot(PCAScore(:,1), PCAScore(:,2),'k.')
% xlabel('PC1')
% ylabel('PC2')

axes(h.p2h)
plot(PCAScore(:,1), PCAScore(:,3),'k.')
% xlabel('PC1')
% ylabel('PC3')
 
axes(h.p3h)
plot(PCAScore(:,1), PCAScore(:,4),'k.')
% xlabel('PC1')
% ylabel('PC4')

axes(h.p4h)
plot(PCAScore(:,2), PCAScore(:,3),'k.')
% xlabel('PC2')
% ylabel('PC3')

axes(h.p5h)
plot(PCAScore(:,2), PCAScore(:,4),'k.')
% xlabel('PC2')
% ylabel('PC4')

axes(h.p6h)
plot(PCAScore(:,3), PCAScore(:,4),'k.')
% xlabel('PC3')
% ylabel('PC4')

set(h.gui, 'Visible','on')

%% Callbacks for guih

    function PC12(hObject, eventdata)
        PCx = 1;
        PCy = 2;
        close;
    end

    function PC13(hObject, eventdata)
        PCx = 1;
        PCy = 3;
        close;
    end

    function PC14(hObject, eventdata)
        PCx = 1;
        PCy = 4;
        close;
    end

    function PC23(hObject, eventdata)
        PCx = 2;
        PCy = 3;
        close;
    end

    function PC24(hObject, eventdata)
        PCx = 2;
        PCy = 4;
        close;
    end

    function PC34(hObject, eventdata)
        PCx = 3;
        PCy = 4;
        close;
    end

    function clusterAll(hObject, eventdata)
        PCx = 0;
        PCy = 0;
        close;
    end

uiwait %waits until the window is closed (as a result of pressing one of the pushbuttons)

end
