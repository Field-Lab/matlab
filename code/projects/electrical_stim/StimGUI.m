function StimGUI

x = 4
%Gui structure
h.fig = figure('position', [400 300 800 600]); %screen size on right-sided display is 1600 by 1200

h.panelOne = uipanel('Title','Electrode Order','FontSize',12,'Position',[.05,.7,.3,.3]); 
h.radioOne = uicontrol('Parent',h.panelOne,'style','radiobutton','FontSize',14,'Position',[10,100,150,40],'String','Random Order'); 

h.panelTwo = uipanel('Title','Waveform Parameters','FontSize',12,'Position',[.5,.1,.45,.9]); 
%h.editOne = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[120,230,30,30])
h.popupOne = uicontrol('Parent',h.panelTwo,'style','popup','FontSize',12,'String',{'1','2','3','4','5','6','7','8'},'Position',[120,385,60,30],'callback',{@addEdit,h}); 
h.textOne = uicontrol('Parent',h.panelTwo,'style','text','FontSize',12,'String','Number of Parts:','Position',[10,390,100,25]); 
h.textTwo = uicontrol('Parent',h.panelTwo,'style','text','FontSize',12,'String','Ratios:','Position',[10,350,100,25]); 
h.editTwo = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[120,350,30,30]); 
%h.editThree = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[120,140,30,30])
%h.editFour = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[120,100,30,30])
h.textThree = uicontrol('Parent',h.panelTwo,'style','text','FontSize',12,'String','Times:','Position',[150,350,100,25]); 
h.editFive = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[250,350,30,30]); 
%h.editSix = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[250,140,30,30])
%h.editSeven = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[250,100,30,30])
h.textFour = uicontrol('Parent',h.panelTwo,'style','text','FontSize',12,'String','Max Value Range:','Position',[10,480,110,25]); 
h.editEight = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[120,480,60,30],'String','Start'); 
h.editNine = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[190,480,60,30],'String','Timestep'); 
h.editTen = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[260,480,60,30],'String','End'); 
h.textFive = uicontrol('Parent',h.panelTwo,'style','text','FontSize',12,'String','Iterations:','Position',[10,430,110,30]); 
h.editEleven = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',12,'Position',[120,430,30,30]); 

h.panelThree = uipanel('Title','Output','FontSize',12,'Position',[.05,.5,.3,.2]); 
h.buttonFive = uicontrol('Parent',h.panelThree,'style','pushbutton','FontSize',12,'String','Choose Save Location','Position',[5,50,150,30]); 
h.editEleven = uicontrol('Parent',h.panelThree,'style','edit','FontSize',12,'String','None Selected','Position',[5,20,150,20],'enable','inactive'); 
disp(h)
disp(guihandles(h.fig))
%Callback functions
function addEdit(hObject, eventdata, h)
	disp
