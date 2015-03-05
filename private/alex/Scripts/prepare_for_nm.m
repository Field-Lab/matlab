%% Prepare Firing rate
clear
cd('/mnt/muench_data/user/alexandra/scripts')

codeWord='NatMov'
dates=cell(1,2);
dates{1}='20130302';
dates{2}='20130302_1';
names=cell(1,39+35);
FR=zeros(27000,144,39+35);
path2save='/mnt/muench_data/data/alexandra/MEA_data/nm12/';
if ~exist(path2save,'dir')
    mkdir(path2save);
end
cnt1=1;
for ttt=dates
    date=cell2mat(ttt)
    
    mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
    heka=dir(['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/*.phys']);
    
    % make list of files to put into easy formatted unit (choosing by the name of the stimulus in heka files)
    file_list=[];
    for i=1:length(heka)
        if ~isempty(regexp(heka(i).name,codeWord, 'once'))
            file_list=[file_list i];
        end
    end
    
    %get list of units
    basic_format_units_list=dir([mainpath, 'units/*.mat']);

    
    for cnt=1:length(basic_format_units_list)

        load([mainpath,'units/',basic_format_units_list(cnt).name]);
        for i=1:length(file_list)
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            a=convolved(spikes,40,27000);
            FR(1:size(a,2)-240,i,cnt1)=a(121:end-120);
            name=basic_format_units_list(cnt).name(1:end-4);
            if cnt1>=40
                name=[name(1:9),'_1',name(10:end)];
            end
            names{cnt1}=name;
        end
        cnt1=cnt1+1;
    end
end
save([path2save,'_FR_NM'],'FR','names')


% plots
figure
p=1;
cnt=1;
clear a
for i=1:6:144

    a(1:27000,p)=mean(FR(:,i:i+5,68),2);
    p=p+1;
    if p==4
        subplot(2,4,cnt)
        plot(a)
        p=1;
        cnt=cnt+1;
    end
    
end



path2save='/mnt/muench_data/data/alexandra/MEA_data/nm12/NM_FR_plot/';
if ~exist(path2save,'dir')
    mkdir(path2save);
end
goodUnits=xlsread('/mnt/muench_data/data/alexandra/MEA_data/nm12/ListOfGoodFiltersPerND.xls');
onOff=goodUnits(:,9);
goodUnits=goodUnits(:,1:8);

figure

for i=1:74
    if ~isempty(regexp(names{i},'0020'))
        i
    end
end
clear a
for kk=1:74
    p=1;
    cnt=1;
    for i=1:6:144
        
        a(1:27000,p)=mean(FR(:,i:i+5,kk),2);
        p=p+1;
        if p==4
            b(1:27000,cnt)=mean(a');
            p=1;
            cnt=cnt+1;
        end
    end
    subplot(1,1,1)
    plot(b(1:26000,2:8))
    spontFR=round((mean(b(1:1500,2:8))*10))/10;
    legend({['ND7  ', num2str(spontFR(1))],['ND6  ', num2str(spontFR(2))],...
        ['ND5  ', num2str(spontFR(3))],['ND4  ', num2str(spontFR(4))],...
        ['ND3  ', num2str(spontFR(5))],['ND2  ', num2str(spontFR(6))],...
        ['ND1  ', num2str(spontFR(7))]})
    mm=max(b(:));
    axis([0 26000 0 mm+1])
    drawnow
    subplot('Position',[0.5 0.96 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    if onOff(kk)>0
        pol='ON';
    elseif onOff(kk)<0
        pol='OFF';
    else
        pol='BAD';
    end
    title([names{kk}, '    ',pol],'FontSize',12,'FontWeight','bold','Interpreter','None')
    drawnow
    saveas(gcf,[path2save,names{kk},'.png']);
end



clear spontFR
for kk=1:74
    p=1;
    cnt=1;
    for i=1:6:144
        
        a(1:27000,p)=mean(FR(:,i:i+5,kk),2);
        p=p+1;
        if p==4
            b(1:27000,cnt)=mean(a');
            p=1;
            cnt=cnt+1;
        end
    end

    spontFR(kk,1:8)=mean(b(1:1500,1:8));

end

figure
subplot(2,1,1)
plot(spontFR(onOff>0,:)')
subplot(2,1,2)
plot(spontFR(onOff<0,:)')


mean(spontFR(onOff>0,:))
std(spontFR(onOff>0,:))
mean(spontFR(onOff<0,:))
std(spontFR(onOff<0,:))


addpath(genpath('/mnt/muench_data/user/alexandra/chronux')) 



figure
clear a
for kk=15
    p=1;
    for i=1:18:144        
        a(1:27000,p)=mean(FR(:,i:i+17,kk),2);
        astd(1:27000,p)=std(FR(:,i:i+17,kk),0,2);
        p=p+1;
    end
    subplot(1,1,1)
    plot(a(:,[3 5]),'linewidth',2)
    hold on
    plot(a(:,[3 5])-astd(:,[3 5]),'--')
    plot(a(:,[3 5])+astd(:,[3 5]),'--')
end


clear a
    clear val lag maxcor maxind
for kk=1:74
    kk
    p=1;
    for i=[37 73]
        a(1:27000,p)=mean(FR(:,i:i+17,kk),2);
        p=p+1;
    end
    cnt=1;

    for i=1:10:26000
        [val(1:201) lag(1:201)]=xcorr(a(i:i+300,1),a(i:i+300,2),100,'coeff');
        [maxcor(kk,cnt),maxind(kk,cnt)]=max(val(:));  
        cnt=cnt+1;
    end
% [val lag]=xcorr(a(:,1),a(:,2),2500,'coeff');
% figure
% 
%     plot(a,'linewidth',2)
%     hold on
%     plot(1:10:26000,lag(maxind)/5,'r.','markersize',20)
%     
%     plot(lag(maxind),'linewidth',2)
%     figure
%     plot(maxcor,'linewidth',2)
%     
%         plot(lag(maxind),maxcor,'*')

end


figure
for kk=1:74
    tmp=lag(maxind(kk,:));
    tmp(tmp==0)=[];
    tmp(tmp==100)=[];
    tmp(tmp==-100)=[];
    [vals(:,kk) bins]=hist([-100 tmp 100],50);
    if onOff(kk)>0
        col='r';
        subplot(3,1,1)
    elseif onOff(kk)<0
        col='b';
        subplot(3,1,2)
    else
        col='k';
        subplot(3,1,3)
    end
    plot(bins,vals(:,kk),col)
    hold on    
end

subplot(3,1,1)
tmp=mean(vals(:,onOff>0),2);
plot(bins,tmp,'r','linewidth',3)
subplot(3,1,2)
tmp=mean(vals(:,onOff<0),2);
plot(bins,tmp,'b','linewidth',3)
subplot(3,1,3)
tmp=mean(vals(:,onOff==0),2);
plot(bins,tmp,'k','linewidth',3)

subplot(3,1,1)
plot(bins,vals(:,15),'m','linewidth',3)
subplot(3,1,2)
plot(bins,vals(:,68),'c','linewidth',3)
        
figure
tmp=mean(vals(:,onOff>0),2);
plot(bins,tmp,'r','linewidth',3)
hold on
tmp=mean(vals(:,onOff<0),2);
plot(bins,tmp,'b','linewidth',3)
tmp=mean(vals(:,onOff==0),2);
plot(bins,tmp,'k','linewidth',3)


plot(lag(maxind(kk,:)),'linewidth',2)
m=maxind(maxind~=101);
hist(lag(m),25)




figure
for kk=1:74
    tmp=maxcor(kk,:);
    [vals(:,kk) bins]=hist([0 tmp 1],50);
    if onOff(kk)>0
        col='r';
        subplot(3,1,1)
    elseif onOff(kk)<0
        col='b';
        subplot(3,1,2)
    else
        col='k';
        subplot(3,1,3)
    end
    plot(bins,vals(:,kk),col)
    hold on    
end
