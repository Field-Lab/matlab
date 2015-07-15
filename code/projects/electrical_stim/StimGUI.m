function StimGUI

%Gui structure
h.fig = figure('position', [400 300 800 600]) %screen size on right-sided display is 1600 by 1200

h.panelOne = uipanel('Title','Electrode Order','FontSize',12,'Position',[.05,.5,.3,.5])
h.radioOne = uicontrol('Parent',h.panelOne,'style','radiobutton','FontSize',14,'Position',[10,180,150,40],'String','Random Order')

h.panelTwo = uipanel('Title','Waveform Parameters','FontSize',12,'Position',[.5,.5,.45,.5])
h.editOne = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[120,230,30,30])
h.textOne = uicontrol('Parent',h.panelTwo,'style','text','FontSize',12,'String','Number of Parts:','Position',[10,230,100,25])
h.textTwo = uicontrol('Parent',h.panelTwo,'style','text','FontSize',12,'String','Ratios:','Position',[10,180,100,25])
h.editTwo = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[120,180,30,30])
h.editThree = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[120,140,30,30])
h.editFour = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[120,100,30,30])
h.textThree = uicontrol('Parent',h.panelTwo,'style','text','FontSize',12,'String','Times:','Position',[150,180,100,25])
h.editFive = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[250,180,30,30])
h.editSix = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[250,140,30,30])
h.editSeven = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[250,100,30,30])
h.textFour = uicontrol('Parent',h.panelTwo,'style','text','FontSize',12,'String','Max Value Range:','Position',[10,50,110,25])
h.editEight = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[120,50,60,30],'String','Start')
h.editNine = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[190,50,60,30],'String','Timestep')
h.editTen = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[260,50,60,30],'String','End')
h.textFive = uicontrol('Parent',h.panelTwo,'style','text','FontSize',12,'String','Iterations:','Position',[10,5,110,30])
h.editEleven = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[120,10,30,30])

%Callback functions
