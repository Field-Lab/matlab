
clear
cd('S:\user\alexandra\scripts')
date='20110925';

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
    if ~isempty(regexp(heka(i).name,'contrast_4','once'))&&isempty(regexp(heka(i).name,'spont_60000', 'once'))
        stim(cnt,1)=20;
        stim(cnt,2)=i;
        cnt=cnt+1;
    elseif ~isempty(regexp(heka(i).name,'contrast_5','once'))&&isempty(regexp(heka(i).name,'spont_60000', 'once'))
        stim(cnt,1)=40;
        stim(cnt,2)=i;
        cnt=cnt+1;
    end
end

%get list of units to be reformatted
basic_format_units_list=dir([mainpath, 'units\*.mat']);
nds='87654322345432'
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
    k=1;
    for i=1:2:12
        subplot(2,6,i)
        plot(mean(all_convolved(cnt_t+1:2:cnt_t+148,:)),'linewidth',2)
        axis([2500,6000,-Inf,Inf])
        title(['ND',nds(k),' RGB 40'])
        subplot(2,6,i+1)
        plot(mean(all_convolved(cnt_t:2:cnt_t+148,:)),'r','linewidth',2)
        axis([2500,6000,-Inf,Inf])
        cnt_t=cnt_t+149;
        title(['ND',nds(k),' RGB 20'])
        k=k+1;
    end 
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title(basic_format_units_list(cnt).name(1:end-3),'Interpreter','none','FontSize',15,'FontWeight','bold')
    saveas(gcf,[path2save,basic_format_units_list(cnt).name(1:end-3),'emf'])
    close(gcf)

end

    