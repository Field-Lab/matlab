function handles=creatFigure(unitName,pol,currentUnit,nUnits,nTrials,start_trial,end_trial,trialsPerScreen)


h=findobj('Name',date);
if ~isempty(h)
    close(h)
end

myFig=figure;
set(myFig,'Name',date,'Units','normalized','Position',[0.0 0.0314 1 0.894])

% change trials per screen
handles.enterTrialsPerScreen=uicontrol('style','edit','String',num2str(trialsPerScreen),'Units','normalized','position',[0.07 0.05 0.02 0.025]);
handles.TrialsPerScreen=uicontrol('style','pushbutton','String','trials/screen','Units','normalized',...
    'position',[0.02 0.05 0.046 0.028],'BackgroundColor',[0.4 0.4 0.9],'Callback','changeTrialsPerScreen(str2num(get(handles.enterTrialsPerScreen,''String'')),currentLFs,start_trial)');

% Polarity Handling
handles.onOff=uibuttongroup('visible','on','Units','normalized','Position',[0.03 0.80 .04 0.08]);
handles.on=uicontrol('style','Radio','String','ON /\','FontWeight','bold','Units','normalized','position',[0.02 0.1 0.9 0.25],...
    'parent',handles.onOff,'HandleVisibility','on','Value',pol);
handles.off=uicontrol('style','Radio','String','OFF \/','FontWeight','bold','Units','normalized','position',[0.02 0.5 0.9 0.25],...
    'parent',handles.onOff,'HandleVisibility','on','Value',~pol);
set(handles.onOff,'SelectionChangeFcn','calcParams(currentLFs,SpikeCount,handles,trialsPerScreen,start_trial)');


% Change center and borders for fitting
handles.changeCenter=uicontrol('style','pushbutton','string','change center');
set(handles.changeCenter,'Units', 'normalized','Position',[0.03 0.75 0.07 0.02],'fontsize',12, 'Callback','change_id(currentLFs,1)');

handles.changeBorders=uicontrol('style','pushbutton','string','change borders');
set(handles.changeBorders,'Units', 'normalized','Position',[0.03 0.7 0.07 0.02],'fontsize',12, 'Callback','change_id(currentLFs,0)');

% make fit/refit
handles.go=uicontrol('style','pushbutton','string','GO!');
set(handles.go,'Units', 'normalized','Position',[0.03 0.65 0.04 0.03],'Backgroundcolor',[1 0.2 0.2],'fontsize',12,...
    'Callback','make_fit(currentLFs, start_trial,trialsPerScreen,fitParams,currentUnit)');

handles.refit=uicontrol('style','pushbutton','string','refit');
set(handles.refit,'Units', 'normalized','Position',[0.03 0.45 0.06 0.03],'Backgroundcolor',[0.6 0.6 0.1],'fontsize',12, 'Callback',...
    'refit(currentLFs, fitParams,currentUnit)');

% go to next trials and units
handles.nextTrials=uicontrol('style','pushbutton','string','next trials');
set(handles.nextTrials,'Units', 'normalized','Position',[0.03 0.60 0.06 0.03],'Backgroundcolor',[0.2 1 0.2],'fontsize',12,...
    'Callback','next_data(0,LinearFilter,SpikeCount,trialsPerScreen,handles,names,polarity)');

handles.nextUnit=uicontrol('style','pushbutton','string','next unit');
set(handles.nextUnit,'Units', 'normalized','Position',[0.03 0.5 0.06 0.03],'Backgroundcolor',[0.1 0.7 0.7],'fontsize',12,...
    'Callback','next_data(1,LinearFilter,SpikeCount,trialsPerScreen,handles,names,polarity)');

handles.enterGotoUnit=uicontrol('style','edit','String',num2str(currentUnit),'Units','normalized','position',[0.07 0.1 0.02 0.025]);
handles.gotoUnit=uicontrol('style','pushbutton','String','go to unit','Units','normalized',...
    'position',[0.02 0.1 0.046 0.028],'BackgroundColor',[0.9 0.6 0.4],...
    'Callback','next_data(2,LinearFilter,SpikeCount,trialsPerScreen,handles,names,polarity)');

% save results
handles.save=uicontrol('style','pushbutton','string','save results');
set(handles.save,'Units', 'normalized','Position',[0.03 0.55 0.08 0.03],'Backgroundcolor',[0.2 0.2 1],'fontsize',12, 'Callback',...
    'save([savePath,''_fitParams''],''fitResults'')');


% print info about the unit (common title)
handles.info=uicontrol('style','text','string',...
    [unitName, '  (',int2str(currentUnit),' out of ',int2str(nUnits),')   Trials ',int2str(start_trial),' to ',int2str(end_trial),' out of ',int2str(nTrials)],...
    'HorizontalAlignment','Left','fontsize',16,'fontweight','bold',...
    'Units', 'normalized','Position',[0.2 0.95 0.7 0.02],'BackgroundColor',[0.8 0.8 0.8]);

% reject one or certain trials
handles.reject=uicontrol('style','pushbutton','string','reject');
set(handles.reject,'Units', 'normalized','Position',[0.03 0.20 0.06 0.03],'Backgroundcolor',[0.9 0.4 0.6],'fontsize',12, 'Callback',...
    'reject_trial(handles,fitParams,currentUnit,currentLFs)');
handles.from=uicontrol('style','edit','String','0','Units','normalized','position',[0.03 0.15 0.02 0.025]);
handles.to=uicontrol('style','edit','String','0','Units','normalized','position',[0.07 0.15 0.02 0.025]);
textFromUI=uicontrol('style','text','String','from','Units','normalized',...
    'position',[0.01 0.15 0.02 0.02],'BackgroundColor',[0.8 0.8 0.8]);
textToUI=uicontrol('style','text','String','to','Units','normalized',...
    'position',[0.051 0.15 0.015 0.02],'BackgroundColor',[0.8 0.8 0.8]);


% set lower fit limits
handles.getFitLimits(1,1)=uicontrol('style','edit','String','0','Units','normalized','position',[0.02 0.35 0.02 0.025]);
handles.getFitLimits(2,1)=uicontrol('style','edit','String','40','Units','normalized','position',[0.02 0.30 0.02 0.025]);
handles.getFitLimits(3,1)=uicontrol('style','edit','String','0','Units','normalized','position',[0.02 0.25 0.02 0.025]);

% set upper fit limits
handles.getFitLimits(1,2)=uicontrol('style','edit','String','1000','Units','normalized','position',[0.08 0.35 0.02 0.025]);
handles.getFitLimits(2,2)=uicontrol('style','edit','String','200','Units','normalized','position',[0.08 0.30 0.02 0.025]);
handles.getFitLimits(3,2)=uicontrol('style','edit','String','100','Units','normalized','position',[0.08 0.25 0.02 0.025]);


textLowerUI=uicontrol('style','text','String','Lower','Units','normalized',...
    'position',[0.02 0.38 0.02 0.025],'BackgroundColor',[0.8 0.8 0.8],'fontweight','bold');
textUpperUI=uicontrol('style','text','String','Upper','Units','normalized',...
    'position',[0.08 0.38 0.02 0.025],'BackgroundColor',[0.8 0.8 0.8],'fontweight','bold');

textAUI=uicontrol('style','text','String','a','Units','normalized',...
    'position',[0.01 0.36 0.01 0.01],'BackgroundColor',[0.8 0.8 0.8],'fontweight','bold');
textBUI=uicontrol('style','text','String','b','Units','normalized',...
    'position',[0.01 0.31 0.01 0.01],'BackgroundColor',[0.8 0.8 0.8],'fontweight','bold');
textCUI=uicontrol('style','text','String','c','Units','normalized',...
    'position',[0.01 0.26 0.01 0.01],'BackgroundColor',[0.8 0.8 0.8],'fontweight','bold');






