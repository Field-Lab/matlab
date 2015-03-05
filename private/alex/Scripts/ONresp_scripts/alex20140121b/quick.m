clear
date='20140121b'

codeWord='quick_';

mainpath=['/Users/alexth/Desktop/old_stuff/',date,'/'];
hekapath=['/Users/alexth/Desktop/old_stuff/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

%find nd changes
nd_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,'nd%nd', 'once'))
        nd_list(1,end+1)=i;        
        a=regexp(heka(i).name,'%nd_');
        a=heka(i).name(a+4);
        nd_list(2,end)=str2num(a);
    end
end

% select heka files
ndFullList=[];
file_list=[];
cnt=1;
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];        
        ndFullList=[ndFullList nd_list(2,find(nd_list(1,:)<file_list(end),1,'last'))];        
    end
end


% get protocol times
protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
protocol=round(protocol(2:end-2,[1,4]));

units=dir([mainpath,'units/*.mat']);

white_flash=zeros(4500,length(file_list),length(units));
black_flash=zeros(4500,length(file_list),length(units));
names=cell(length(units),1);
for cnt=1:length(units)
    
    load([mainpath,'units/',units(cnt).name]);
    
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        conv=convolved(spikes,40,75000);
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash(:,i,cnt)=white_flash(:,i,cnt)+conv(protocol(st,1)-499:protocol(st,1)+4000)'/5;
            black_flash(:,i,cnt)=black_flash(:,i,cnt)+conv(protocol(st+2,1)-499:protocol(st+2,1)+4000)'/5;
        end
    end  
    
    if isempty(regexp(units(cnt).name,date, 'once'))
        names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
    else
        names{cnt}=units(cnt).name(1:end-4);
    end
end

save([path2save,'quick'],'ndFullList','white_flash','black_flash','names')


%% Plot

clear

date='20140121b'
path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
load([path2load,'quick'])

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/quick_plot/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

nds='87654';
figure
set(gcf,'position',[388    87   557   854])

for j=1:size(black_flash,3)
    
    for cnt=1:5
        subplot(5,1,cnt)
        hold off
        plot(black_flash(:,ndFullList==str2num(nds(cnt)),j))
        hold on
        plot(4701:9200,white_flash(:,ndFullList==str2num(nds(cnt)),j))
        title(['ND',nds(cnt)])
        a=get(gca,'YLim');
        line([500, 500],[0 a(2)],'color','k')
        line([2500, 2500],[0 a(2)],'color','k')
        line([5200, 5200],[0 a(2)],'color','k')
        line([7200, 7200],[0 a(2)],'color','k')
    end
  
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{j},'  Black, White'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{j},'.png'])
end

