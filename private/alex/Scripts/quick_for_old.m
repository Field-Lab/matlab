clear

cd('/mnt/muench_data/user/alexandra/scripts')
date='20110925'
codeWord_black='contrast_4'; % black
codeWord='contrast_5'; % black

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

% select heka files
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))&isempty(regexp(heka(i).name,'spont_60000', 'once'))&isempty(regexp(heka(i).name,'sh_2', 'once'))
        file_list=[file_list i];
    end
end
file_list_black=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord_black, 'once'))&isempty(regexp(heka(i).name,'spont_60000', 'once'))&isempty(regexp(heka(i).name,'sh_2', 'once'))
        file_list_black=[file_list_black i];
    end
end

% get protocol times
protocol=read_header_field_heka(hekapath, heka(file_list(50)).name, 'Stimulus Protocol');
protocol=round(protocol(1:end,[1,2]));
beg=2951;

units=dir([mainpath,'units/*.mat']);
white_flash=zeros(6000,length(file_list),length(units));
black_flash=zeros(6000,length(file_list),length(units));
names=cell(length(units),1);
for cnt=1:length(units)
    
    load([mainpath,'units/',units(cnt).name]);
    
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        conv=convolved(spikes,40,12500);
        conv=conv(121:end-120);
        white_flash(:,i,cnt)=conv(beg-1999:beg+4000)';
    end
    
    for i=1:length(file_list_black)
        spikes=round(cell2mat(unit{1,2}(file_list_black(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        conv=convolved(spikes,40,12500);
        conv=conv(121:end-120);
        black_flash(:,i,cnt)=conv(beg-1999:beg+4000)';
    end
    
    if isempty(regexp(units(cnt).name,date, 'once'))
        names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
    else
        names{cnt}=units(cnt).name(1:end-4);
    end
end

save([path2save,'quick'],'white_flash','black_flash','names')

wf=zeros(6000,14,length(units));
bf=zeros(6000,14,length(units));
cc=1;
for i=1:74:1036
    wf(:,cc,:)=reshape(mean(white_flash(:,i+54:i+73,:),2),6000,47);
    bf(:,cc,:)=reshape(mean(black_flash(:,i+54:i+73,:),2),6000,47);
    cc=cc+1;
end
save([path2save,'quick'],'wf','bf','names')

%% Plot every cell separately
path2save='/mnt/muench_data/data/alexandra/MEA_data/analysis/20110925/QUICK_Pictures/';
if ~exist(path2save,'dir')
    mkdir(path2save)
end
figure
nds='87654322'
for cnt=1:47
    maxFR=[];
    for i=1:8
        subplot(4,2,i)
        hold off
        plot(wf(:,i,cnt),'r','linewidth',2) 
        hold on
        plot(6201:12200, bf(:,i,cnt),'linewidth',2)        
        maxFR=[maxFR max(wf(:,i,cnt)) max(bf(:,i,cnt))];      
    end
    k=max(maxFR)*1.1;
    for i=1:8
        subplot(4,2,i)
        axis([0 12200 0 k])
        line([2000,2000],[0,k],'color','k')
        line([4000,4000],[0,k],'color','k')
        line([8200,8200],[0,k],'color','k')
        line([10200,10200],[0,k],'color','k')
        title(['ND',nds(i)])
    end
    subplot('Position',[0.5 0.96 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{cnt},'   RED - WHITE FLASH, BLUE - BLACK FLASH'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{cnt},'.png'])    
end


