function StimGUI
	%Got help from: http://blogs.mathworks.com/videos/2007/12/28/matlab-basics-guis-without-guide/

	%define background variables
	clear all
	global numParts
    global editFields
	global ratioFields
	global timeFields
	global h
	numParts = 1;
	editFields = {};
	ratioFields = {};
	timeFields = {};

	%get screen size
	set(0,'units','pixels');
	ss = get(0,'screensize');
	sw = ss(3);
	sh = ss(4);

	%set fonts
	header_font = 12;
	reg_font = 8; 

	%Gui structure
	h.fig = figure('position', [.3*sw, .15*sh, .5*sw, .65*sh]); %screen size on right-sided display is 1600 by 1200

	%Panel One--------------------
	h.panelOne = uipanel('Title','Electrode Order','FontSize',header_font,'Position',[.05,.7,.45,.3]); 
	h.p1radioOne = uicontrol('Parent',h.panelOne,'style','radiobutton','FontSize',reg_font,'Units','normalized','Position',[.1,.6,.6,.1],'String','Random Order'); 
	h.p1checkOne = uicontrol('Parent',h.panelOne,'style','checkbox','FontSize',reg_font,'Units','normalized','Position',[.6,.6,.4,.1],'String','Random Amplitude Order');

	%Panel Two--------------------
	h.panelTwo = uipanel('Title','Waveform Parameters','FontSize',header_font,'Position',[.5,.1,.45,.9]); 
	%Panel Parameters (must be put here because for some reason if defined after particular element value is not accessible)
	h.p2textheight = .04;
	h.p2textwidth = .3;
	h.p2g1voffset = .92; 
	h.p2g1editw = .15; 
	h.p2g1editheight = .05; 
	h.p2g2voffset = .82;
	h.p2g2popupheight = .05; 
	h.p2g2popupwidth = .15; 
	h.p2g3voffset = .7;
	h.p2g3editwidth = .1;
	h.p2g3editheight = .05;
	h.p2g3releditvoffset = .1
	%first row 
	h.p2textOne = uicontrol('Parent',h.panelTwo,'style','text','FontSize',reg_font,'String','Max Value Range:','Units','normalized','Position',[0,h.p2g1voffset,h.p2textwidth,h.p2textheight]); 
	h.p2g1editOne = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',reg_font,'Units','normalized','Position',[.35,h.p2g1voffset,h.p2g1editw,h.p2g1editheight],'String','Start'); 
	h.p2g1editTwo = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',reg_font,'Units','normalized','Position',[.55,h.p2g1voffset,h.p2g1editw,h.p2g1editheight],'String','Timestep'); 
	h.p2g1editThree = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',reg_font,'Units','normalized','Position',[.75,h.p2g1voffset,h.p2g1editw,h.p2g1editheight],'String','End'); 
	%second row
	h.p2g2textOne = uicontrol('Parent',h.panelTwo,'style','text','FontSize',reg_font,'String','Number of Parts:','Units','normalized','Position',[0,h.p2g2voffset,h.p2textwidth,h.p2textheight]); 
	%h.p2g2popupOne = uicontrol('Parent',h.panelTwo,'style','popup','FontSize',reg_font,'String',{'1','2','3','4','5','6','7','8'},'Units','normalized','Position',[.35,h.p2g2voffset,h.p2g2popupwidth,h.p2g2popupheight],'callback',{@addEdit,h}); 
	h.p2g2popupOne = uicontrol('Parent',h.panelTwo,'style','popup','FontSize',reg_font,'String',{'1','2','3','4','5','6','7','8'},'Units','normalized','Position',[.35,h.p2g2voffset,h.p2g2popupwidth,h.p2g2popupheight],'callback',@addEdit); 
	h.p2g2textTwo = uicontrol('Parent',h.panelTwo,'style','text','FontSize',reg_font,'String','Iterations:','Units','normalized','Position',[.5,h.p2g2voffset,h.p2textwidth,h.p2textheight]);
	h.p2g2editOne = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',reg_font,'Units','normalized','Position',[.75,h.p2g2voffset,h.p2g3editwidth,h.p2g3editheight]);
	%remaining rows
	h.p2g3textOne = uicontrol('Parent',h.panelTwo,'style','text','FontSize',reg_font,'String','Ratios:','Units','normalized','Position',[0,h.p2g3voffset,h.p2textwidth,h.p2textheight]); 
	h.p2g3editOne = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',reg_font,'Units','normalized','Position',[.35,h.p2g3voffset,h.p2g3editwidth,h.p2g3editheight]); 
	h.p2g3textTwo = uicontrol('Parent',h.panelTwo,'style','text','FontSize',reg_font,'String','Times:','Units','normalized','Position',[.5,h.p2g3voffset,h.p2textwidth,h.p2textheight]); 
	h.p2g3editTwo = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',reg_font,'Units','normalized','Position',[.75,h.p2g3voffset,h.p2g3editwidth,h.p2g3editheight]); 

	%Panel Three--------------------
	h.panelThree = uipanel('Title','Output','FontSize',header_font,'Position',[.05,.3,.45,.4]); 
	h.p3buttonOne = uicontrol('Parent',h.panelThree,'style','pushbutton','FontSize',reg_font,'String','Choose Save Location','Units','normalized','Position',[.01,.8,.6,.15],'callback',@chooseSave); 
	h.p3editOne = uicontrol('Parent',h.panelThree,'style','edit','FontSize',reg_font,'String','None Selected','Units','normalized','Position',[.01,.675,.6,.1]); 
	h.p3buttonTwo = uicontrol('Parent',h.panelThree,'style','pushbutton','FontSize',reg_font+2,'String','Go!','Units','normalized','Position',[.7,.675,.2,.2],'BackgroundColor',[.2344,.6992,.4414],'callback',@goButton); 
	h.p3textOne = uicontrol('Parent',h.panelThree,'style','text','FontSize',reg_font,'String','Filename:','Units','normalized','Position',[0,.4,.35,.1])
	h.p3editTwo = uicontrol('Parent',h.panelThree,'style','edit','FontSize',reg_font,'String','jklfs','Units','normalized','Position',[.4,.4,.35,.15]); 

	%Callback functions
	function addEdit(hObject, eventdata) %ISSUE: Currently deletes contents of cell when chaging number
		prev_numParts = numParts;
		numParts = get(hObject, 'Value');
		%Remove previous fields
		if (prev_numParts ~= numParts) && (~isempty(editFields))
			for i=editFields
				delete(h.(i{1}));
				h = rmfield(h,i);
			end
			editFields = {};
			ratioFields = {};
			timeFields = {};
		end
		%add new fields
		if numParts > 1
			for i=2:numParts
				newEditField = strcat('p2g3editratio', num2str(i));
				editFields = [editFields newEditField];
				ratioFields = [ratioFields newEditField]; %just collect ratio field handle
				h.(newEditField) = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',reg_font,'Units','normalized','Position',[.35,h.p2g3voffset-h.p2g3releditvoffset*(i-1),h.p2g3editwidth,h.p2g3editheight]);
				newEditField = strcat('p2g3edittime', num2str(i));
				editFields = [editFields newEditField];
				timeFields = [timeFields newEditField]; %just collect time field handle
				h.(newEditField) = uicontrol('Parent',h.panelTwo,'style','edit','FontSize',reg_font,'Units','normalized','Position',[.75,h.p2g3voffset-h.p2g3releditvoffset*(i-1),h.p2g3editwidth,h.p2g3editheight]);
			end
		end
	end

	function chooseSave(hObject,evendata)
		dn = uigetdir;
		set(h.p3editOne,'String',dn);
	end
	function goButton(hObject,eventdata)
		%Checks and Warnings by order of uipanel--------------------
		%Electrode and Amplitude Order
		randflg = get(h.p1radioOne,'Value');
		if randflg ~= 1; errordlg('Please Specify Electrode Order'); return; end
		randampflg = 0; if (get(h.p1checkOne,'Value') == get(h.p1checkOne,'Max')); randampflg = 1; end %get checkbox value
		%Current Value Range and Step
		currBegin = str2num(get(h.p2g1editOne,'String'));
		currStep = str2num(get(h.p2g1editTwo,'String'));
		currEnd = str2num(get(h.p2g1editThree,'String'));
		if (isempty(currBegin)) || (isempty(currStep)) || (isempty(currEnd))
			errordlg('One of the Max Value Range entries is empty or non-numeric');
			return;
		elseif (currEnd - currBegin) < currStep
			errordlg('Timestep is too large');
			return;
		end
		%Iterations
		iter = str2num(get(h.p2g2editOne,'String'));
		if (isempty(iter)) || (iter == 0); errordlg('Iterations is non-numeric or 0'); return; end
		%Ratios and Times (numparts needn't be collected because it should be up to date)
		%get first part
		r = []; t = [];
		r = [r str2num(get(h.p2g3editOne,'String'))];
		t = [t str2num(get(h.p2g3editTwo,'String'))];
		if (isempty(r)) || (isempty(t)); errordlg('Must specify at least one part'); return; end
		%iterate through rest
		for i=ratioFields; r = [r str2num(get(h.(i{1}),'String'))]; end
		for i=timeFields; t = [t str2num(get(h.(i{1}),'String'))]; end
		for i=t; if rem(i/50,1) ~= 0; errordlg('Time increments must be multiples of 50 us'); return ; end; end %make sure all time increments are multiples of 50
		%Check save location
		dn = get(h.p3editOne,'String'); 
		if strcmp(dn,'None Selected') || strcmp(dn,'0'); errordlg('Must select save location'); return; end
		fn = get(h.p3editTwo,'String'); %get file base name
		%Call Processing Function (should be in same directory as GUI)
		gapp = 7500;
		processGuiOutput([currBegin currStep currEnd], iter, r, t, randflg,randampflg,gapp,fn,dn);
		msgbox('Done!');
	end
end
