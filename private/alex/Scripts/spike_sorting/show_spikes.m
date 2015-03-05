function show_spikes(flag)
% shows units mean and std. flag: 1 - show overlay of selected units,
% 0 - show separately in subplots

global rgc spikes uicontr_units
persistent sbpl_cnt h
col=[1,0,0;0,1,0;0.5,0.5,0;0.25,0.75,0.75;0.75,0.25,0.75;0.05,0.58,0.05;0.9,0.45,0.1;0.4,0.25,0.75;0.86,0.86,0];
try
    if flag
        unit2overlay=[];
        for i=1:size(uicontr_units,1)
            if get(uicontr_units(i,1),'Value')==1
                unit2overlay=[unit2overlay i];
            end
        end
        if isempty(unit2overlay)
            unit2overlay=1:size(uicontr_units,1);
        end
        figure('Units','normalized','position',[0.002 0.0343 0.2 0.4075],'Name','Units Overlay')
        set(subplot(1,1,1),'Color','k')
        hold on
        for i=1:length(unit2overlay)
            try
                mean_Unit=mean(spikes(:,rgc{unit2overlay(i)}),2);
                std_Unit=std(spikes(:,rgc{unit2overlay(i)}),0,2);
            catch err
                display('ERROR! (memory?) Split spikes. STD is not to be trusted')
                mean_Unit=mean(spikes(:,rgc{unit2overlay(i)}(1:round(end/2))),2);
                mean_Unit=mean_Unit+mean(spikes(:,rgc{unit2overlay(i)}(round(end/2)+1:end)),2);
            end
            patch([1:24 24:-1:1],[mean_Unit+std_Unit; mean_Unit(end:-1:1)-std_Unit(end:-1:1)],col(mod(unit2overlay(i)-1,9)+1,:),'FaceAlpha',0.7,'EdgeAlpha',0.7);
            plot(1:24,mean_Unit,'color',col(mod(unit2overlay(i)-1,9)+1,:),'LineWidth',2);
        end
        title(['UNITS ',int2str(unit2overlay)])
        hold off
    else
        h = findobj('Name','Units mean');
        if ~isempty(h)
            for i=1:length(sbpl_cnt)
                delete(sbpl_cnt(i))
            end
            clear sbpl_cnt
        else
            h=figure('Units','normalized','position',[0.206 0.0343 0.307 0.4075],'Name','Units mean');
        end
        figure(h)
        % amount of units (max 9)
        if size(rgc,2)<=9
            unit2show=size(rgc,2);
        else
            display('You really, really wanna more than 9 units?... code a new figure ^^')
            unit2show=9;
        end
        % subplot grid
        n=ceil(unit2show/3);
        if unit2show==1
            l=1;
        elseif unit2show==2 || unit2show==4
            l=2;
        else
            l=3;
        end
        
        for k=1:unit2show
            sbpl_cnt(k)=subplot(n,l,k);
            set(sbpl_cnt(k),'Color','k');
            hold on
            try
                mean_Unit=mean(spikes(:,rgc{k}),2);
                std_Unit=std(spikes(:,rgc{k}),0,2);
            catch err
                display('ERROR! (memory?) Split spikes. STD is not to be trusted')
                mean_Unit=mean(spikes(:,rgc{k}(1:round(end/2))),2);
                mean_Unit=mean_Unit+mean(spikes(:,rgc{k}(round(end/2)+1:end)),2);
            end
            min_Unit=min(spikes(:,rgc{k})');
            max_Unit=max(spikes(:,rgc{k})');
            patch([1:24 24:-1:1],[mean_Unit+std_Unit; mean_Unit(end:-1:1)-std_Unit(end:-1:1)],col(mod(k-1,9)+1,:),'FaceAlpha',0.7,'EdgeAlpha',0.7);
            plot(1:24,mean_Unit,'color',col(mod(k-1,9)+1,:),'LineWidth',2);
            plot(1:24,min_Unit,'color',col(mod(k-1,9)+1,:),'LineStyle','--','LineWidth',1);
            plot(1:24,max_Unit,'color',col(mod(k-1,9)+1,:),'LineStyle','--','LineWidth',1);
            title(['Unit',int2str(k)],'FontSize',10)
            axis tight
        end
    end
catch err
    display(err)
end
