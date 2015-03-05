
clear
cd('S:\user\alexandra\scripts')
date='20110809';

mainpath=['S:\data\alexandra\MEA_data\',date,'\'];
hekapath=['S:\data\alexandra\MEA_data\',date,'\HEKA\'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath 'pictures\'],'dir')    
    mkdir([mainpath, 'pictures']);
end
path2save=[mainpath 'pictures\'];

% prepare stimulus
cnt=1;
for i=1:length(heka)-2
    if ~isempty(regexp(heka(i).name,'contrast_2','once'))&&isempty(regexp(heka(i).name,'spont_60000', 'once'))
        stim(cnt,1)=120;
        stim(cnt,2)=i;
        cnt=cnt+1;
    elseif ~isempty(regexp(heka(i).name,'contrast_3','once'))&&isempty(regexp(heka(i).name,'spont_60000', 'once'))
        stim(cnt,1)=215;
        stim(cnt,2)=i;
        cnt=cnt+1;
    end
end

%get list of units to be reformatted
basic_format_units_list=dir([mainpath, 'units\*.mat']);

for cnt=1:length(basic_format_units_list)
    basic_format_units_list(cnt).name
    load([mainpath,'units\',basic_format_units_list(cnt).name]);
    
    all_convolved=zeros(length(stim),8000);
    for i=1:length(stim)
        tmp=round(cell2mat(unit{1,2}(stim(i,2),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        tmp_conv=convolved(tmp,40,8000);
        all_convolved(i,1:8000)=tmp_conv(121:end-120);
    end    
    
    set(gcf,'Position', [1 31 1680 946])
    cnt_t=1;
    for i=1:9
        subplot(3,3,i)
        plot(mean(all_convolved(cnt_t+1:2:cnt_t+148,:)),'linewidth',2)
        hold on
        plot(mean(all_convolved(cnt_t:2:cnt_t+148,:)),'r','linewidth',2)
        cnt_t=cnt_t+149;
        title(['ND',int2str(10-i)])
    end 
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title(basic_format_units_list(cnt).name(1:end-3),'Interpreter','none','FontSize',15,'FontWeight','bold')
    saveas(gcf,[path2save,basic_format_units_list(cnt).name(1:end-3),'emf'])
    close(gcf)

end

  set(gcf,'Position', [1 31 1680 946])
    cnt_t=1;
    for i=1:2
        subplot(3,3,i)
        plot(mean(all_convolved(cnt_t+1:2:cnt_t+148,:)),'linewidth',2)
        hold on
        plot(mean(all_convolved(cnt_t:2:cnt_t+148,:)),'r','linewidth',2)
        cnt_t=cnt_t+149;
        title(['ND',int2str(10-i)])
    end 
    
    k=1
    for i=cnt_t+1:10:cnt_t+148
        subplot(5,3,k)
        plot(mean(all_convolved(i:2:i+8,:)))
        k=k+1;        
    end
