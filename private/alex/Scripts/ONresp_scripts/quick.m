clear
date='2013-09-03'
date='20130906'
% date='20130723a'
date='20130703a'
date='20130703b'

codeWord='quick_';

mainpath=['/Users/alexth/Desktop/old_stuff/',date,'/'];
hekapath=['/Users/alexth/Desktop/old_stuff/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

% select heka files
clear ndFullList
file_list=[];
cnt=1;
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
        
        a=regexp(heka(i).name,'ND');
        a=heka(i).name(a+2);
        if length(a)>1
            ndFullList(cnt)=str2num(a(1))+str2num(a(2));
        else
            ndFullList(cnt)=str2num(a(1));
        end
        cnt=cnt+1;
    end
end


% get protocol times
protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
protocol=round(protocol(2:end-2,[1,2]));
protocol(protocol(1:end,2)==1,2)=50;
protocol(protocol(1:end,2)==0,2)=10;
protocol(protocol(1:end,2)==-1,2)=30;

units=dir([mainpath,'units/*.mat']);

% if strcmp(date,'20130723a')
%     for cnt=1:length(units)
%         load([mainpath,'units/',units(cnt).name]);
%         for i=1:length(heka)
%             unit{1,2}{i,1}=heka(i).name;            
%         end
%         save([mainpath,'units/',units(cnt).name],'unit');
%     end
% end

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


%% Plot for Ringer Switch

clear

date='20130723a'
path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
load([path2load,'quick'])

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/quick_plot/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end


for j=1:size(black_flash,3)
    figure
    set(gcf,'position',[388    87   557   854])
    
    subplot(5,1,5)
    plot(black_flash(:,1:8,j),'b')
    hold on
    plot(4701:9200,white_flash(:,1:8,j),'b')
    title('ND8')
    
    subplot(5,1,4)
    plot(black_flash(:,9:12,j),'b')
    hold on
    plot(black_flash(:,13:16,j),'r')
    plot(4701:9200,white_flash(:,9:12,j),'b')
    plot(4701:9200,white_flash(:,13:16,j),'r')
    title('ND7')
    
    subplot(5,1,3)
    plot(black_flash(:,17:22,j),'b')
    hold on
    plot(black_flash(:,23:26,j),'r')
    plot(4701:9200,white_flash(:,17:22,j),'b')
    plot(4701:9200,white_flash(:,23:26,j),'r')
    title('ND6')
    
    subplot(5,1,2)
    plot(black_flash(:,27:30,j),'b')
    hold on
    plot(black_flash(:,31:34,j),'r')
    plot(4701:9200,white_flash(:,27:30,j),'b')
    plot(4701:9200,white_flash(:,31:34,j),'r')
    title('ND5')
    
    subplot(5,1,1)
    plot(black_flash(:,35:40,j),'b')
    hold on
    plot(black_flash(:,41:44,j),'r')
    plot(4701:9200,white_flash(:,35:40,j),'b')
    plot(4701:9200,white_flash(:,41:44,j),'r')
    title('ND4')
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{j},'  Black, White'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{j},'.png'])
    close all
end




%% Plot for APB

clear

date='2013-09-03'
path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
load([path2load,'quick'])

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/quick_plot/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end


for j=1:size(black_flash,3)
    figure
    set(gcf,'position',[4          87        1243         849])
    
    subplot(2,2,1)
    plot(black_flash(:,1:6,j),'b')
    hold on
    plot(4701:9200,white_flash(:,1:6,j),'b')
    title('ND7')
    
    subplot(2,2,2)
    plot(black_flash(:,7:10,j),'b')
    hold on
    plot(black_flash(:,11:18,j),'r')
    plot(black_flash(:,19:30,j),'g')
    plot(4701:9200,white_flash(:,7:10,j),'b')
    plot(4701:9200,white_flash(:,11:18,j),'r')
    plot(4701:9200,white_flash(:,19:30,j),'g')
    title('ND6')
    
    subplot(2,2,3)
    plot(black_flash(:,31:34,j),'b')
    hold on
    plot(black_flash(:,35:42,j),'r')
    plot(black_flash(:,43:54,j),'g')
    plot(4701:9200,white_flash(:,31:34,j),'b')
    plot(4701:9200,white_flash(:,35:42,j),'r')
    plot(4701:9200,white_flash(:,43:54,j),'g')
    title('ND5')
    
    subplot(2,2,4)
    plot(black_flash(:,55:58,j),'b')
    hold on
    plot(black_flash(:,59:66,j),'r')
    plot(black_flash(:,67:74,j),'g')
    plot(4701:9200,white_flash(:,55:58,j),'b')
    plot(4701:9200,white_flash(:,59:66,j),'r')
    plot(4701:9200,white_flash(:,67:74,j),'g')
    title('ND4')
    
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{j},'  Black White'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{j},'.png'])
    close all
end





figure
set(gcf,'position',[4          87        1243         849])

for j=1:size(black_flash,3)
    
    subplot(2,2,1)
    hold off
    plot(mean(black_flash(:,3:6,j),2),'b','linewidth',2)
    hold on
    plot(4701:9200,mean(white_flash(:,3:6,j),2),'b','linewidth',2)
    title('ND7')
    
    subplot(2,2,2)
    hold off
    plot(mean(black_flash(:,7:10,j),2),'b','linewidth',2)
    hold on
    plot(mean(black_flash(:,15:18,j),2),'r','linewidth',2)
    plot(mean(black_flash(:,27:30,j),2),'g','linewidth',2)
    plot(4701:9200,mean(white_flash(:,7:10,j),2),'b','linewidth',2)
    plot(4701:9200,mean(white_flash(:,15:18,j),2),'r','linewidth',2)
    plot(4701:9200,mean(white_flash(:,27:30,j),2),'g','linewidth',2)
    title('ND6')
    
    subplot(2,2,3)
    hold off
    plot(mean(black_flash(:,31:34,j),2),'b','linewidth',2)
    hold on
    plot(mean(black_flash(:,39:42,j),2),'r','linewidth',2)
    plot(mean(black_flash(:,51:54,j),2),'g','linewidth',2)
    plot(4701:9200,mean(white_flash(:,31:34,j),2),'b','linewidth',2)
    plot(4701:9200,mean(white_flash(:,39:42,j),2),'r','linewidth',2)
    plot(4701:9200,mean(white_flash(:,51:54,j),2),'g','linewidth',2)
    title('ND5')
    
    subplot(2,2,4)
    hold off
    plot(mean(black_flash(:,55:58,j),2),'b','linewidth',2)
    hold on
    plot(mean(black_flash(:,63:66,j),2),'r','linewidth',2)
    plot(mean(black_flash(:,75:78,j),2),'g','linewidth',2)
    legend('C','APB','OUT')
    plot(4701:9200,mean(white_flash(:,55:58,j),2),'b','linewidth',2)
    plot(4701:9200,mean(white_flash(:,63:66,j),2),'r','linewidth',2)
    plot(4701:9200,mean(white_flash(:,75:78,j),2),'g','linewidth',2)
    title('ND4')
    
    
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{j},'  Black White'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{j},'.png'])
end
close all



% special for 20130723a

figure
set(gcf,'position',[4          87        1243         849])

for j=1:size(black_flash,3)
    
    subplot(2,2,1)
    hold off
    plot(mean(black_flash(:,3:6,j),2),'b','linewidth',2)
    hold on
    plot(4701:9200,mean(white_flash(:,3:6,j),2),'b','linewidth',2)
    title('ND7')
    
    subplot(2,2,2)
    hold off
    plot(mean(black_flash(:,7:10,j),2),'b','linewidth',2)
    hold on
    plot(mean(black_flash(:,15:18,j),2),'r','linewidth',2)
    plot(mean(black_flash(:,24:27,j),2),'g','linewidth',2)
    plot(4701:9200,mean(white_flash(:,7:10,j),2),'b','linewidth',2)
    plot(4701:9200,mean(white_flash(:,15:18,j),2),'r','linewidth',2)
    plot(4701:9200,mean(white_flash(:,24:27,j),2),'g','linewidth',2)
    title('ND6')
    
    subplot(2,2,3)
    hold off
    plot(mean(black_flash(:,28:31,j),2),'b','linewidth',2)
    hold on
    plot(mean(black_flash(:,36:39,j),2),'r','linewidth',2)
    plot(mean(black_flash(:,48:51,j),2),'g','linewidth',2)
    plot(4701:9200,mean(white_flash(:,28:31,j),2),'b','linewidth',2)
    plot(4701:9200,mean(white_flash(:,36:39,j),2),'r','linewidth',2)
    plot(4701:9200,mean(white_flash(:,48:51,j),2),'g','linewidth',2)
    title('ND5')
    
    subplot(2,2,4)
    hold off
    plot(mean(black_flash(:,52:55,j),2),'b','linewidth',2)
    hold on
    plot(mean(black_flash(:,60:63,j),2),'r','linewidth',2)
    plot(mean(black_flash(:,72:75,j),2),'g','linewidth',2)
    legend('C','APB','OUT')
    plot(4701:9200,mean(white_flash(:,52:55,j),2),'b','linewidth',2)
    plot(4701:9200,mean(white_flash(:,60:63,j),2),'r','linewidth',2)
    plot(4701:9200,mean(white_flash(:,72:75,j),2),'g','linewidth',2)
    title('ND4')
    
    
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{j},'  Black White'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{j},'.png'])
end
close all




% special for 2013-09-03

figure
set(gcf,'position',[4          87        1243         849])

for j=1:size(black_flash,3)
    
    subplot(2,2,1)
    hold off
    plot(mean(black_flash(:,3:6,j),2),'b','linewidth',2)
    hold on
    plot(4701:9200,mean(white_flash(:,3:6,j),2),'b','linewidth',2)
    title('ND7')
    
    subplot(2,2,2)
    hold off
    plot(mean(black_flash(:,7:10,j),2),'b','linewidth',2)
    hold on
    plot(mean(black_flash(:,15:18,j),2),'r','linewidth',2)
    plot(mean(black_flash(:,27:30,j),2),'g','linewidth',2)
    plot(4701:9200,mean(white_flash(:,7:10,j),2),'b','linewidth',2)
    plot(4701:9200,mean(white_flash(:,15:18,j),2),'r','linewidth',2)
    plot(4701:9200,mean(white_flash(:,27:30,j),2),'g','linewidth',2)
    title('ND6')
    
    subplot(2,2,3)
    hold off
    plot(mean(black_flash(:,31:34,j),2),'b','linewidth',2)
    hold on
    plot(mean(black_flash(:,39:42,j),2),'r','linewidth',2)
    plot(mean(black_flash(:,51:54,j),2),'g','linewidth',2)
    plot(4701:9200,mean(white_flash(:,31:34,j),2),'b','linewidth',2)
    plot(4701:9200,mean(white_flash(:,39:42,j),2),'r','linewidth',2)
    plot(4701:9200,mean(white_flash(:,51:54,j),2),'g','linewidth',2)
    title('ND5')
    
    subplot(2,2,4)
    hold off
    plot(mean(black_flash(:,55:58,j),2),'b','linewidth',2)
    hold on
    plot(mean(black_flash(:,63:66,j),2),'r','linewidth',2)
    plot(mean(black_flash(:,77:80,j),2),'g','linewidth',2)
    legend('C','APB','OUT')
    plot(4701:9200,mean(white_flash(:,55:58,j),2),'b','linewidth',2)
    plot(4701:9200,mean(white_flash(:,63:66,j),2),'r','linewidth',2)
    plot(4701:9200,mean(white_flash(:,77:80,j),2),'g','linewidth',2)
    title('ND4')
    
    
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{j},'  Black White'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{j},'.png'])
end
close all
