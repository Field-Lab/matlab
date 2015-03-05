function create_figure(filtersFile,flag)

global fitLines fitParams 

if flag==1 % first call
    
    load(filtersFile)
    
    h=findobj('Name',date);
    if ~isempty(h)
        close(h)
    end
    
    tmp=strfind(filtersFile,filesep);
    savePath=filtersFile(1:tmp(end));
    
    trialsPerScreen=min(size(LinearFilter,2),49);
    currentUnit=1;
    currentPolarity=polarity(currentUnit);
    pol=(currentPolarity/2)+0.5;
    start_trial=1;
    end_trial=trialsPerScreen;
    
    mainGUI=figure;
    set(mainGUI,'Name',date,'Units','normalized','Position',[0.0 0.0314 1 0.894])
    
    onOffUI = uibuttongroup('visible','on','Units','normalized','Position',[0.03 0.80 .04 0.08]);
    onUI=uicontrol('style','Radio','String','ON /\','FontWeight','bold','Units','normalized','position',[0.02 0.1 0.9 0.25],...
        'parent',onOffUI,'HandleVisibility','on','Value',pol);
    offUI=uicontrol('style','Radio','String','OFF \/','FontWeight','bold','Units','normalized','position',[0.02 0.5 0.9 0.25],...
        'parent',onOffUI,'HandleVisibility','on','Value',~pol);
    set(onOffUI,'SelectionChangeFcn','pol_change(currentLFs,get(onUI,''Value''),sbpl,start_trial)');
    
    changeIndUI=uicontrol('style','pushbutton','string','change index');
    set(changeIndUI,'Units', 'normalized','Position',[0.03 0.75 0.07 0.02],'fontsize',12, 'Callback','change_id(1)');
    
    changeIntUI=uicontrol('style','pushbutton','string','change borders');
    set(changeIntUI,'Units', 'normalized','Position',[0.03 0.7 0.07 0.02],'fontsize',12, 'Callback','change_id(2)');
    
    goUI=uicontrol('style','pushbutton','string','GO!');
    set(goUI,'Units', 'normalized','Position',[0.03 0.65 0.04 0.03],'Backgroundcolor',[1 0.2 0.2],'fontsize',12, 'Callback','make_fit(toFit, sbpl, fitParams, getFitLimits, start_trial, end_trial)');
    
    nextUI=uicontrol('style','pushbutton','string','next trials');
    set(nextUI,'Units', 'normalized','Position',[0.03 0.60 0.06 0.03],'Backgroundcolor',[0.2 1 0.2],'fontsize',12, 'Callback','fittingGUI(0)');
    
    saveUI=uicontrol('style','pushbutton','string','save results');
    set(saveUI,'Units', 'normalized','Position',[0.03 0.55 0.08 0.03],'Backgroundcolor',[0.2 0.2 1],'fontsize',12, 'Callback',...
        'save([path2save,unit_name,''_fit_res_HC.mat''],''common_res_fit'')');
    
    nextUnitUI=uicontrol('style','pushbutton','string','next unit');
    set(nextUnitUI,'Units', 'normalized','Position',[0.03 0.5 0.06 0.03],'Backgroundcolor',[0.1 0.7 0.7],'fontsize',12, 'Callback',...
        'fittingGUI(1)');
    
    infoUI=uicontrol('style','text','string',...
        [names{currentUnit}, '  (',int2str(currentUnit),' out of ',int2str(size(LinearFilter,3)),')   Trials ',int2str(start_trial),' to ',int2str(end_trial),' out of ',int2str(size(LinearFilter,2))],...
        'HorizontalAlignment','Left','fontsize',16,'fontweight','bold',...
        'Units', 'normalized','Position',[0.2 0.95 0.7 0.02],'BackgroundColor',[0.8 0.8 0.8]);
    
    refitUI=uicontrol('style','pushbutton','string','refit');
    set(refitUI,'Units', 'normalized','Position',[0.03 0.45 0.06 0.03],'Backgroundcolor',[0.6 0.6 0.1],'fontsize',12, 'Callback',...
        'refit()');
    
    % lower fit limits
    getFitLimits(1,1)=uicontrol('style','edit','String','0','Units','normalized','position',[0.02 0.35 0.02 0.025]);
    getFitLimits(2,1)=uicontrol('style','edit','String','40','Units','normalized','position',[0.02 0.30 0.02 0.025]);
    getFitLimits(3,1)=uicontrol('style','edit','String','0','Units','normalized','position',[0.02 0.25 0.02 0.025]);
    
    %upper fit limits
    getFitLimits(1,2)=uicontrol('style','edit','String','1000','Units','normalized','position',[0.08 0.35 0.02 0.025]);
    getFitLimits(2,2)=uicontrol('style','edit','String','200','Units','normalized','position',[0.08 0.30 0.02 0.025]);
    getFitLimits(3,2)=uicontrol('style','edit','String','100','Units','normalized','position',[0.08 0.25 0.02 0.025]);
    
    
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
    
    rejectUI=uicontrol('style','pushbutton','string','reject');
    set(rejectUI,'Units', 'normalized','Position',[0.03 0.20 0.06 0.03],'Backgroundcolor',[0.9 0.4 0.6],'fontsize',12, 'Callback',...
        'reject_trial()');
    fromUI=uicontrol('style','edit','String','0','Units','normalized','position',[0.03 0.15 0.02 0.025]);
    toUI=uicontrol('style','edit','String','0','Units','normalized','position',[0.07 0.15 0.02 0.025]);
    textFromUI=uicontrol('style','text','String','from','Units','normalized',...
        'position',[0.01 0.15 0.02 0.02],'BackgroundColor',[0.8 0.8 0.8]);
    textToUI=uicontrol('style','text','String','to','Units','normalized',...
        'position',[0.051 0.15 0.015 0.02],'BackgroundColor',[0.8 0.8 0.8]);
    
    gotoUnitUI=uicontrol('style','edit','String',num2str(currentUnit),'Units','normalized','position',[0.07 0.1 0.02 0.025]);
    textGotoUnitUI=uicontrol('style','pushbutton','String','go to unit','Units','normalized',...
        'position',[0.02 0.1 0.046 0.028],'BackgroundColor',[0.9 0.6 0.4],'Callback','fittingGUI(2)');
    
    trialsPerScreenUI=uicontrol('style','edit','String',num2str(trialsPerScreen),'Units','normalized','position',[0.07 0.05 0.02 0.025]);
    textGotoUnitUI=uicontrol('style','pushbutton','String','trials/screen','Units','normalized',...
        'position',[0.02 0.05 0.046 0.028],'BackgroundColor',[0.4 0.4 0.9],'Callback','fittingGUI(3)'); 
end

if start_trial==1
    
    currentLFs=LinearFilter(:,:,currentUnit);
    
    if pol %on cell
        [amp,ind]=max(currentLFs(51:200,:));
    else
        [amp,ind]=min(currentLFs(51:200,:));
    end
    fitParams.ind=ind+50;
    fitParams.rightBorder=fitParams.ind+50;
    fitParams.leftBorder=fitParams.ind-50;
    fitParams.pol=pol;
    fitParams.amp=amp;
    fitParams.currentPolarity=currentPolarity;
    fitParams.SpikeCount=SpikeCount(:,currentUnit);
end

[rows,cols]=opt_subplots(trialsPerScreen);
toFit=currentLFs(:,start_trial:end_trial);
sbpl=zeros(end_trial-start_trial+1,1);
for i=1:end_trial-start_trial+1
    trial=start_trial+i-1;
    sbpl(i)=subplot(rows,cols,i);
    hold off
    plot(toFit(:,i));
    axis tight    
    fitLines(2,i)=line([fitParams.ind(trial),fitParams.ind(trial)],[min(toFit(:,i)),max(toFit(:,i))],'color','k','linewidth',2);
    fitLines(1,i)=line([fitParams.rightBorder(trial),fitParams.rightBorder(trial)],[min(toFit(:,i)),max(toFit(:,i))],'color','g','linewidth',2);
    fitLines(3,i)=line([fitParams.leftBorder(trial),fitParams.leftBorder(trial)],[min(toFit(:,i)),max(toFit(:,i))],'color','g','linewidth',2);
    title(['trial ', int2str(trial),' (',int2str(i),'), ',int2str(SpikeCount(trial,currentUnit)),' sp']);
end

