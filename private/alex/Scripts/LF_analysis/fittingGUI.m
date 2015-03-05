function fittingGUI(filtersFile,flag)


create_figure(flag)

    


% fitting GUI
global HighFilters HighSpikeCount LowFilters LowSpikeCount spike_info mainGUI...
    toFit sbpl line_id ind rightBorder leftBorder rightBorder_id leftBorder_id...    
    start_trial cur_unit pol fl path2save unit_name getLowerAUI getIniAUI...
    getUpperAUI getLowerBUI getIniBUI getUpperBUI getLowerCUI getIniCUI getUpperCUI...
    onUI offUI common_res_fit spikeCount fromUI toUI...
    filterspath filters unitspath units goto date gotoUnitUI ...
    trialsPerScreen trialsPerScreenUI infoUI onOffUI changeIndUI changeIntUI ...
    goUI nextUI saveUI nextUnitUI refitUI rejectUI AllFilters AllTrials



    
    unitspath=[pathName,'easy_formatted_units/'];

    cur_unit=1;
    start_trial=1;    
    fl=1;
    if strcmp(contr,'high')
        filterspath=[pathName,'FFFlicker_LF/'];
        filters=dir([filterspath,'*filter.mat']);
        load([filterspath,filters(1).name])
        units=dir([unitspath,'*FFFlicker.mat']);
        path2save=[pathName,'fit_res_HC/'];
        AllTrials=size(HighFilters,2);
    elseif strcmp(contr,'low')
        filterspath=[pathName,'FFFlicker_LF/'];
        filters=dir([filterspath,'*filter.mat']);
        load([filterspath,filters(1).name])
        units=dir([unitspath,'*FFFlicker.mat']);
        path2save=[pathName,'fit_res_LC/'];
        AllTrials=size(LowFilters,2);
    elseif strcmp(contr,'checker')
        filterspath=[pathName,'CBFlicker_LF/'];
        filters=dir([filterspath,'*filter.mat']);
        load([filterspath,filters(1).name])
        units=dir([unitspath,'*CBFlicker.mat']);
        path2save=[pathName,'fit_res_CB/'];
        AllTrials=size(HighFilters,2);
    end
    if ~exist(path2save,'dir')
        mkdir(path2save);
    end
    
    h=findobj('Name',date);
    if ~isempty(h)
        close(h)
    end
    mainGUI=figure;
    set(mainGUI,'Name',date,'Units','normalized','Position',[0.0 0.0314 1 0.894])
    
    onOffUI = uibuttongroup('visible','on','Units','normalized','Position',[0.03 0.80 .04 0.08]);
    onUI=uicontrol('style','Radio','String','ON /\','FontWeight','bold','Units','normalized','position',[0.02 0.1 0.9 0.25],...
        'parent',onOffUI,'HandleVisibility','on','Value',ceil(pol+0.5)/2);
    offUI=uicontrol('style','Radio','String','OFF \/','FontWeight','bold','Units','normalized','position',[0.02 0.5 0.9 0.25],...
        'parent',onOffUI,'HandleVisibility','on','Value',~(ceil(pol+0.5)/2));
    set(onOffUI,'SelectionChangeFcn','pol_change()');
    
    changeIndUI=uicontrol('style','pushbutton','string','change index');
    set(changeIndUI,'Units', 'normalized','Position',[0.03 0.75 0.07 0.02],'fontsize',12, 'Callback','change_id(1)');
    
    changeIntUI=uicontrol('style','pushbutton','string','change borders');
    set(changeIntUI,'Units', 'normalized','Position',[0.03 0.7 0.07 0.02],'fontsize',12, 'Callback','change_id(2)');
    
    goUI=uicontrol('style','pushbutton','string','GO!');
    set(goUI,'Units', 'normalized','Position',[0.03 0.65 0.04 0.03],'Backgroundcolor',[1 0.2 0.2],'fontsize',12, 'Callback','make_fit()');
    
    nextUI=uicontrol('style','pushbutton','string','next trials');
    set(nextUI,'Units', 'normalized','Position',[0.03 0.60 0.06 0.03],'Backgroundcolor',[0.2 1 0.2],'fontsize',12, 'Callback','fittingGUI(0)');
 
    saveUI=uicontrol('style','pushbutton','string','save results');
    set(saveUI,'Units', 'normalized','Position',[0.03 0.55 0.08 0.03],'Backgroundcolor',[0.2 0.2 1],'fontsize',12, 'Callback',...
        'save([path2save,unit_name,''_fit_res_HC.mat''],''common_res_fit'')');
    
    nextUnitUI=uicontrol('style','pushbutton','string','next unit');
    set(nextUnitUI,'Units', 'normalized','Position',[0.03 0.5 0.06 0.03],'Backgroundcolor',[0.1 0.7 0.7],'fontsize',12, 'Callback',...
        'fittingGUI(1)');
    
    infoUI=uicontrol('style','text','string',...
        ['Unit:  ', unit_name, '  (',int2str(cur_unit),' out of ',int2str(length(filters)),')    Trials:  ',int2str(start_trial),'-',int2str(goto),'  (out of ',int2str(AllTrials),')    ',contr,' Contrast'],...
        'HorizontalAlignment','Left','fontsize',16,'fontweight','bold',...
        'Units', 'normalized','Position',[0.2 0.95 0.7 0.02],'BackgroundColor',[0.8 0.8 0.8]);
    
    refitUI=uicontrol('style','pushbutton','string','refit');
    set(refitUI,'Units', 'normalized','Position',[0.03 0.45 0.06 0.03],'Backgroundcolor',[0.6 0.6 0.1],'fontsize',12, 'Callback',...
        'refit()');
    
    getLowerAUI=uicontrol('style','edit','String','0','Units','normalized','position',[0.02 0.35 0.02 0.025]);
    getLowerBUI=uicontrol('style','edit','String','40','Units','normalized','position',[0.02 0.30 0.02 0.025]);
    getLowerCUI=uicontrol('style','edit','String','0','Units','normalized','position',[0.02 0.25 0.02 0.025]);
    
    getIniAUI=uicontrol('style','edit','String','100','Units','normalized','position',[0.05 0.35 0.02 0.025]);
    getIniBUI=uicontrol('style','edit','String','100','Units','normalized','position',[0.05 0.30 0.02 0.025]);
    getIniCUI=uicontrol('style','edit','String','30','Units','normalized','position',[0.05 0.25 0.02 0.025]);
    
    getUpperAUI=uicontrol('style','edit','String','1000','Units','normalized','position',[0.08 0.35 0.02 0.025]);
    getUpperBUI=uicontrol('style','edit','String','200','Units','normalized','position',[0.08 0.30 0.02 0.025]);
    getUpperCUI=uicontrol('style','edit','String','100','Units','normalized','position',[0.08 0.25 0.02 0.025]);

    
    textLowerUI=uicontrol('style','text','String','Lower','Units','normalized',...
        'position',[0.02 0.38 0.02 0.025],'BackgroundColor',[0.8 0.8 0.8],'fontweight','bold');
    textIniUI=uicontrol('style','text','String','Initial','Units','normalized',...
        'position',[0.05 0.38 0.02 0.025],'BackgroundColor',[0.8 0.8 0.8],'fontweight','bold');
    textUpperUI=uicontrol('style','text','String','Upper','Units','normalized',...
        'position',[0.08 0.38 0.02 0.025],'BackgroundColor',[0.8 0.8 0.8],'fontweight','bold');
    
    textAUI=uicontrol('style','text','String','a','Units','normalized',...
        'position',[0.01 0.36 0.01 0.01],'BackgroundColor',[0.8 0.8 0.8],'fontweight','bold');
    textBUI=uicontrol('style','text','String','b','Units','normalized',...
        'position',[0.01 0.31 0.01 0.01],'BackgroundColor',[0.8 0.8 0.8],'fontweight','bold');
    textCUI=uicontrol('style','text','String','c','Units','normalized',...
        'position',[0.01 0.26 0.01 0.01],'BackgroundColor',[0.8 0.8 0.8],'fontweight','bold');
    
    rejectUI=uicontrol('style','pushbutton','string','reject');
    set(rejectUI,'Units', 'normalized','Position',[0.03 0.20 0.06 0.03],'Backgroundcolor',[0.9 0.4 0.6],'fontsize',12, 'Callback',...
        'reject_trial()');
    fromUI=uicontrol('style','edit','String','0','Units','normalized','position',[0.03 0.15 0.02 0.025]);
    toUI=uicontrol('style','edit','String','0','Units','normalized','position',[0.07 0.15 0.02 0.025]);
    textFromUI=uicontrol('style','text','String','from','Units','normalized',...
        'position',[0.01 0.15 0.02 0.02],'BackgroundColor',[0.8 0.8 0.8]);
    textToUI=uicontrol('style','text','String','to','Units','normalized',...
        'position',[0.051 0.15 0.015 0.02],'BackgroundColor',[0.8 0.8 0.8]);
    
    gotoUnitUI=uicontrol('style','edit','String',num2str(cur_unit),'Units','normalized','position',[0.07 0.1 0.02 0.025]);
    textGotoUnitUI=uicontrol('style','pushbutton','String','go to unit','Units','normalized',...
        'position',[0.02 0.1 0.046 0.028],'BackgroundColor',[0.9 0.6 0.4],'Callback','fittingGUI(2)');
    
    trialsPerScreenUI=uicontrol('style','edit','String',num2str(trialsPerScreen),'Units','normalized','position',[0.07 0.05 0.02 0.025]);
    textGotoUnitUI=uicontrol('style','pushbutton','String','trials/screen','Units','normalized',...
        'position',[0.02 0.05 0.046 0.028],'BackgroundColor',[0.4 0.4 0.9],'Callback','fittingGUI(3)');
    
    
elseif ifNextUnit==1
    cur_unit=cur_unit+1;
    start_trial=1;
    fl=1;
elseif ifNextUnit==2    
    cur_unit=str2num(get(gotoUnitUI,'String'));
    start_trial=1;
    fl=1;
elseif ifNextUnit==3
    trialsPerScreen=str2num(get(trialsPerScreenUI,'String'));
else    
    if goto<AllTrials
        start_trial=goto+1;
        fl=fl+1;
    else
        cur_unit=cur_unit+1;
        start_trial=1;
        fl=1;
    end
end
    
if start_trial==1
    load([filterspath,filters(cur_unit).name])
    if strcmp(contr,'high') | strcmp(contr,'checker')
        AllFilters=HighFilters(1:350,:);
        AllTrials=size(HighFilters,2);  
        spikeCount=HighSpikeCount;
    else
        AllFilters=LowFilters(1:350,:);
        AllTrials=size(LowFilters,2);
        spikeCount=LowSpikeCount;
    end
    common_res_fit=zeros(3,AllTrials);
    load([unitspath,units(cur_unit).name])
    unit_name=units(cur_unit).name(1:end-27);
end

goto=min(AllTrials,trialsPerScreen*fl); %last trial to display

if start_trial==1
    myHope=std(AllFilters(50:150,:))-std(AllFilters(1:50,:))>5;
    [~, ind_max]=max(AllFilters(1:200,:));
    [~, ind_min]=min(AllFilters(1:200,:));
    if sum(ind_max(myHope)<ind_min(myHope))<sum(myHope)/2
        pol=-1; % off
        ind=ind_min;
        ind(ind<51)=51;
    else
        pol=1; % on
        ind=ind_max;
        ind(ind<51)=51;
    end    
    rightBorder=ind+50;
    leftBorder=ind-50;
end

set(onUI, 'Value',ceil(pol+0.5)/2) ;
set(offUI,'Value',~ceil(pol+0.5)/2);
set(infoUI,'String', ['Unit:  ', unit_name, '  (',int2str(cur_unit),' out of ',int2str(length(filters)),')    Trials:  ',int2str(start_trial),'-',int2str(goto),'  (out of ',int2str(AllTrials),')    ',contr,' Contrast']);
set(fromUI,'String','0');
set(toUI,'String','0');
set(gotoUnitUI,'String',num2str(cur_unit));
set(trialsPerScreenUI,'String',num2str(trialsPerScreen));


for i=1:length(sbpl)
    if ishandle(sbpl(i))
        delete(sbpl(i));
    end
end

[rows,cols]=opt_subplots(goto-start_trial+1);
toFit=AllFilters(:,start_trial:goto);
sbpl=zeros(goto-start_trial+1,1);
for i=1:goto-start_trial+1
    j=start_trial+i-1;
    sbpl(i)=subplot(rows,cols,i);
    hold off
    plot(toFit(:,i));
    axis tight    
    line_id(i)=line([ind(j),ind(j)],[min(toFit(:,i)),max(toFit(:,i))],'color','k','linewidth',2);
    rightBorder_id(i)=line([rightBorder(j),rightBorder(j)],[min(toFit(:,i)),max(toFit(:,i))],'color','g','linewidth',2);
    leftBorder_id(i)=line([leftBorder(j),leftBorder(j)],[min(toFit(:,i)),max(toFit(:,i))],'color','g','linewidth',2);
    title([int2str(i),' trial ',int2str(j),', ',int2str(spikeCount(j)),' sp']);
end

end