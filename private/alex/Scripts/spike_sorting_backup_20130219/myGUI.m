function myGUI
global spikes pc1 pc2 tdata rgc spikeTimes isi_all uicontr_units src hekaNames traceUI ...
    getSlice1 getSlice2 thr spikes2delete getXdata getYdata otherSignUI getFatUI...
    qualityUI notesUI autoPathUI unitUI fileName savePath timeUnit keepAxesUI undoUI...
    selectValidUI getTemplateX getTemplateY
    
% default: reset view to original
ifKeepAxes=0;

%initialize buffer
undo(2,1,1);

% pca
x=spikes(:,min(spikes)<thr);
x=cov(x');
[V,~]=eig(x);
pc_vectors=V(:,end-2:end);
pc1=spikes'*pc_vectors(:,end);
pc2=spikes'*pc_vectors(:,end-1);

%time vector
tdata=(1:size(spikes,2))';

%vector for ISI
isi_all=zeros(1,size(spikes,2)+1);
isi_all(1)=0;
k=1;
for i=1:length(spikeTimes)
    isi_all(k+1:k+length(spikeTimes{i}))=spikeTimes{i}+isi_all(k)+1;
    k=k+length(spikeTimes{i});
end
isi_all(1)=[];


%preview figure
h=findobj('Name','Preview');
if ~isempty(h)
    close(h)
end
prv=figure;
set(gcf,'Name','Preview','Units','normalized','Position',[ 0.205    0.07    0.77    0.8])
subplot('Position',[0.03 0.05 0.2 0.9])
plot(pc1(1:10:end),tdata(1:10:end),'.','MarkerSize',0.1)
set(gca,'YTick',0,'yticklabel','')
xlabel('pc1'); ylabel('time (spike N)')
axis tight
subplot('Position',[0.25 0.05 0.2 0.9])
plot(pc2(1:10:end),tdata(1:10:end),'.','MarkerSize',0.1)
set(gca,'YTick',0,'yticklabel','')
xlabel('pc2'); ylabel('time (spike N)')
axis tight
subplot('Position',[0.47 0.05 0.2 0.9])
plot(min(spikes(:,1:10:end)),tdata(1:10:end),'.','MarkerSize',0.1)
set(gca,'YTick',0,'yticklabel','')
xlabel('min'); ylabel('time (spike N)')
axis tight
subplot('Position',[0.71 0.53 0.28 0.42])
plot(min(spikes(:,1:10:end)),max(spikes(:,1:10:end)),'.','MarkerSize',0.1)
xlabel('min'); ylabel('max')
axis tight
subplot('Position',[0.71 0.05 0.28 0.42])
plot(pc1(1:10:end),pc2(1:10:end),'.','MarkerSize',0.1)
xlabel('pc1'); ylabel('pc2')
axis tight


rgc=cell(1,1);
uicontr_units=[];
if ~isempty(src)
    close(src)
end
% h=findobj('Name','Units tracing GUI');
% if ~isempty(h)
%     close(h)
% end
h=findobj('Name','Units mean');
if ~isempty(h)
    close(h)
end
h=findobj('Name','ISI <5ms; % - less than 2ms');
if ~isempty(h)
    close(h)
end

src=figure;
set(gcf,'Name',['Units tracing GUI, ',fileName],'Units','normalized','Position',[0.5167 0.0314 0.4815 0.8943])
ah = axes('DrawMode','fast');
set(gca,'Position',[0.08 0.07 0.6 0.85])
plot(pc1,tdata,'.','MarkerSize',0.1)
xlabel('pc1')
ylabel('time (spike N)')

%set spike features to display
getXdata=uicontrol('style','popup','Units','normalized','position',[0.87 0.89 0.08 0.025],...
    'string', {'pc1' 'pc2' 'slice 1' 'slice 2' 'time' 'min' 'max' 'templ X'},'fontweight','bold', ...
    'callback', 'redraw(get(getXdata,''value''),get(getYdata,''value''))');
set(getXdata,'Value',1)
textXdata=uicontrol('Style','text','BackgroundColor',[0.8 0.8 0.8],'Units','normalized',...
    'position',[0.951 0.885 0.05 0.025],'String','X axis','FontWeight','bold');

getYdata=uicontrol('style','popup','Units','normalized','position',[0.87 0.85 0.08 0.025],...
    'string',{'pc1' 'pc2' 'slice 1' 'slice 2' 'time' 'min' 'max' 'templ Y'},'fontweight','bold',...
    'callback', 'redraw(get(getXdata,''value''),get(getYdata,''value''))');
set(getYdata,'Value',5)
textYdata=uicontrol('Style','text','BackgroundColor',[0.8 0.8 0.8],'Units','normalized',...
    'position',[0.951 0.845 0.05 0.025],'String','Y axis','FontWeight','bold');

getSlice1=uicontrol('style','edit','String','13','Units','normalized','position',[0.7 0.89 0.03 0.025],...    
    'callback', 'redraw(get(getXdata,''value''),get(getYdata,''value''))');
textSlice1=uicontrol('Style','text','BackgroundColor',[0.8 0.8 0.8],'Units','normalized',...
    'position',[0.7305 0.895 0.05 0.015],'String','Slice 1');

getSlice2=uicontrol('style','edit','String','18','Units','normalized','position',[0.7 0.85 0.03 0.025],...   
    'callback', 'redraw(get(getXdata,''value''),get(getYdata,''value''))');
textSlice2=uicontrol('Style','text','BackgroundColor',[0.8 0.8 0.8],'Units','normalized',...
    'Position',[0.7305 0.855 0.05 0.015],'String','Slice 2');

getTemplateX=uicontrol('style','edit','String','1','Units','normalized','position',[0.78 0.89 0.03 0.025],...    
    'callback', 'redraw(get(getXdata,''value''),get(getYdata,''value''))');
textTemplateX=uicontrol('Style','text','BackgroundColor',[0.8 0.8 0.8],'Units','normalized',...
    'position',[0.8105 0.895 0.05 0.015],'String','Templ X');

getTemplateY=uicontrol('style','edit','String','2','Units','normalized','position',[0.78 0.85 0.03 0.025],...   
    'callback', 'redraw(get(getXdata,''value''),get(getYdata,''value''))');
textTemplateY=uicontrol('Style','text','BackgroundColor',[0.8 0.8 0.8],'Units','normalized',...
    'Position',[0.8105 0.855 0.05 0.015],'String','Templ Y');



% manipulate units

% delete unit
deleteUI=uicontrol('style','pushbutton','string','delete units');
set(deleteUI,'Units','normalized','position',[0.7 0.8 0.15 0.03],...
    'fontsize',10,'Callback','manipulate_units(1,get(getXdata,''value''),get(getYdata,''value''))');

% join unit
joinUI=uicontrol('style','pushbutton','string','join units');
set(joinUI,'Units','normalized','position',[0.7 0.76 0.15 0.03],...
    'fontsize',10,'Callback','manipulate_units(3,get(getXdata,''value''),get(getYdata,''value''))');

% make spikes exlusive (belong only to the selected unit)
exclusiveSpikesUI=uicontrol('style','pushbutton','string','make spikes exclusive');
set(exclusiveSpikesUI,'Units','normalized','position',[0.7 0.72 0.2 0.03],...
    'fontsize',10,'Callback','manipulate_units(2,get(getXdata,''value''),get(getYdata,''value''))');

% get unit from overlapping spikes
overlapSpikesUI=uicontrol('style','pushbutton','string','unit from overlap (n=2)');
set(overlapSpikesUI,'Units','normalized','position',[0.7 0.68 0.2 0.03],...
    'fontsize',10,'Callback','manipulate_units(8,get(getXdata,''value''),get(getYdata,''value''))');

% delete spikes
deleteSpikesUI=uicontrol('style','pushbutton','string','delete spikes');
set(deleteSpikesUI,'Units','normalized','position',[0.7 0.63 0.15 0.03],...
    'fontsize',10,'Callback','manipulate_units(4,get(getXdata,''value''),get(getYdata,''value''))');
deleteSpikesTextUI=uicontrol('Style','text','BackgroundColor',[0.8 0.8 0.8],'Units','normalized',...
    'position',[0.86 0.62 0.05 0.05],'String',{'that are','also','in unit'});
spikes2delete=uicontrol('style','edit','Units','normalized','position',[0.93 0.632 0.03 0.025]);

% get template
getTemplateUI=uicontrol('style','pushbutton','string','get template');
set(getTemplateUI,'Units','normalized','position',[0.87 0.8 0.1 0.03],...
    'fontsize',10,'Callback','manipulate_units(9,get(getXdata,''value''),get(getYdata,''value''))');

% template projection
selectValidUI=uicontrol('style','togglebutton','string','select valid','BackgroundColor',[0.2 0.9 0.4],'fontsize',10);
set(selectValidUI,'Units','normalized','position',[0.87 0.76 0.1 0.03],'Callback','manipulate_units(10,get(getXdata,''value''),get(getYdata,''value''))');


% templateUI=uicontrol('style','pushbutton','string','plot template');
% set(templateUI,'Units','normalized','position',[0.87 0.76 0.1 0.03],...
%     'fontsize',10,'Callback','manipulate_units(10,get(getXdata,''value''),get(getYdata,''value''))');

% make invisible
deleteUI=uicontrol('style','pushbutton','string','make (in)visible');
set(deleteUI,'Units','normalized','position',[0.7 0.58 0.15 0.03],...
    'fontsize',10,'Callback','manipulate_units(5,get(getXdata,''value''),get(getYdata,''value''))');

% fat points
fatUI=uicontrol('style','pushbutton','string','make fat');
set(fatUI,'Units','normalized','position',[0.7 0.54 0.15 0.03],...
    'fontsize',10,'Callback','manipulate_units(6,get(getXdata,''value''),get(getYdata,''value''))');
getFatUI=uicontrol('style','edit','Units','normalized','position',[0.87 0.542 0.03 0.025]);

% other sign
signUI=uicontrol('style','pushbutton','string','change sign');
set(signUI,'Units','normalized','position',[0.7 0.5 0.15 0.03],...
    'fontsize',10,'Callback','manipulate_units(7,get(getXdata,''value''),get(getYdata,''value''))');
otherSignUI=uicontrol('style','edit','Units','normalized','position',[0.87 0.502 0.03 0.025]);

% show overlay of units
overlayUI=uicontrol('style','pushbutton','string','units overlay');
set(overlayUI,'Units','normalized','position',[0.7 0.46 0.15 0.03],...
    'fontsize',10,'Callback', 'show_spikes(1)');

% show features maps of selected units
mapUI=uicontrol('style','pushbutton','string','units map');
set(mapUI,'Units','normalized','position',[0.88 0.46 0.10 0.03],...
    'fontsize',10,'Callback', 'show_map');



%basic ui controls: update, forel, trace, get, save, undo
updateUI=uicontrol('style','pushbutton','string','update');
set(updateUI,'Units','normalized','position',[0.87 0.96 0.099 0.03],'fontsize',12, 'Callback','redraw(get(getXdata,''value''),get(getYdata,''value''))');

forelUI=uicontrol('style','pushbutton','string','load forel');
set(forelUI,'Units','normalized','position',[0.81 0.92 0.091 0.025],'fontsize',12, 'Callback','load_forel()');

doForelUI=uicontrol('style','pushbutton','string','do forel');
set(doForelUI,'Units','normalized','position',[0.91 0.92 0.083 0.025],'fontsize',12, 'Callback','do_forel()');

undoUI=uicontrol('style','pushbutton','string','undo','BackgroundColor',[0.6 0.6 0.9],'fontsize',12);
set(undoUI,'Units','normalized','position',[0.7 0.965 0.1 0.02],'fontsize',12,'Callback','undo(1,get(getXdata,''value''),get(getYdata,''value''))');


keepAxesUI=uicontrol('style','togglebutton','string','keep axes','BackgroundColor',[0.9 0.6 0.6],'fontsize',12);
set(keepAxesUI,'Units','normalized','position',[0.7 0.935 0.1 0.02],'fontsize',12);


traceUI=uicontrol('style','togglebutton','String','press to trace',...
    'BackgroundColor',[0.5 0.5 0.5],'fontsize',12,'Units','normalized',...
    'position',[0.725 0.42 0.2 0.03],'Callback','tracing(src,traceUI)');
traceInstrUI=uicontrol('style','text','string',...
    {'Start/end: click within plot','','s - show path / pause tracing (resume if ''c'' pressed afterwards)','','c - resume tracing'},...
    'HorizontalAlignment','Left','fontsize',10,'fontweight','bold',...
    'Units','normalized','position',[0.72 0.1 0.25 0.3],'BackgroundColor',[0.8 0.8 0.8]);

textUnitUI=uicontrol('style','text','string',{'Time Unit (1 - s, 1000 - ms)'},...
    'Units','normalized','position',[0.755 0.21 0.18 0.02],'BackgroundColor',[0.8 0.8 0.8]);
unitUI=uicontrol('style','edit','String',timeUnit,'Units','normalized','position',[0.71 0.21 0.05 0.025]);

textNoteUI=uicontrol('style','text','string',{'notes (free text)'},...
    'Units','normalized','position',[0.7 0.17 0.25 0.03],'BackgroundColor',[0.8 0.8 0.8]);
notesUI=uicontrol('style','edit','Units','normalized','position',[0.71 0.13 0.24 0.05]);

textQualityUI=uicontrol('style','text','string',{'insert quality for all units to save','separation: space, order: increasing'},...
    'Units','normalized','position',[0.7 0.085 0.25 0.03],'BackgroundColor',[0.8 0.8 0.8]);
qualityUI=uicontrol('style','edit','Units','normalized','position',[0.71 0.05 0.24 0.025],...
    'Callback','quality=get(qualityUI,''string'');');

saveUnitsUI=uicontrol('style','pushbutton','string','save units');
set(saveUnitsUI,'Units','normalized','position',[0.75 0.01 0.15 0.03],'fontsize',12,'Callback',...
    'save_units(date,savePath,spikeTimes,thr,thrs,ch,quality,experimenter,hekaNames,[get(getXdata,''value''),get(getYdata,''value'')])');

autoPathUI=uicontrol('style','checkbox','value',1);
set(autoPathUI,'Units','normalized','position',[0.71 0.02 0.02 0.016],'BackgroundColor',[0.8 0.8 0.8]);
autoPathTextUI=uicontrol('style','text','string','Auto Path');
set(autoPathTextUI,'Units','normalized','position',[0.64 0.018 0.07 0.017],'BackgroundColor',[0.8 0.8 0.8]);
autoPathLocationUI=uicontrol('style','text','string',['save to: ',savePath,'units']);
set(autoPathLocationUI,'Units','normalized','position',[0.14 0.002 0.5 0.017],'fontsize',10,'BackgroundColor',[0.8 0.8 0.8]);


set(gcf,'Toolbar','figure')

