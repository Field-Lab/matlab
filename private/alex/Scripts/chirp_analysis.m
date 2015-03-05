date='20121023'

unitsPath=dir(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/easy_formatted_units/*chirp*'])
trials=size(spike_info.name_info,1);
[rows,cols]=opt_subplots(trials*2);
nds='765432123456';
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946]);
for unit=1:length(unitsPath)
    load(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/easy_formatted_units/',unitsPath(unit).name])
    chirp=zeros(trials,25000);
    for trial=1:trials
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end,1));
        conv=conv(121:end-120);
        chirp(trial,1:length(conv))=conv;
    end

    cnt=1;
    
    for nd=6:5:40
        subplot(7,1,cnt)
        hold off
        tmp=mean(chirp(nd:nd+4,1:20000));
        plot((1:length(tmp))/10000,tmp,'LineWidth',1);
        k=mean(tmp(500:1500));
        line([0,length(tmp)/10000],[k,k],'color','m')
        title(['ND',nds(cnt)],'FontSize',12,'FontWeight','bold')
        for i=1:length(c)
            line([c(i),c(i)],[0,max(tmp)],'color','k')
        end        
        hold on
        plot((1:length(tmp))/10000,tmp,'LineWidth',1);
        axis tight
        set(gca,'xtick',0)
        cnt=cnt+1;
    end
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([unitsPath(unit).name(1:end-23),'. Chirp spike rate, average of 5 trials, end of ND'],'FontSize',20,'FontWeight','bold','Interpreter','None')
    saveas(gcf,['/mnt/muench_data/user/alexandra/MEA_data/',date,'/chirp_plots/',unitsPath(unit).name(1:end-23),'_chirp.jpg']);

end

cnt=1;
for nd=1:5:40
    tmp(cnt,1:25000)=mean(chirp(nd:nd+4,:));
    cnt=cnt+1;
end

figure
plot(tmp','linewidth',2)
legend(['87654321']')

tmp2=tmp(:,3330:11200)';
figure
plot(tmp2,'linewidth',2)
legend(['87654321']')

figure
plot(spike_info.flip_times{1}(:,2))
tmp3=spike_info.flip_times{1}(2:end-1,2);
[a,b]=findpeaks(tmp3);
figure;
plot(tmp3)
hold on
plot(b,a,'r*')

c=spike_info.flip_times{1}(b+1,1)
c=c/10000;

tmp2=tmp(:,3300:11200)';
figure
plot(tmp2(:,2:8),'linewidth',2)
for i=1:length(c)
    line([c(i)-3300,c(i)-3300],[0,max(tmp2(:))],'color','k')
end
legend(['7654321']')
f=c-3300;
for i=1:length(f)-1
    for m=1:8
        [a,b]=findpeaks(tmp2(f(i):f(i+1),m));
        if isempty(b)
            b=-1;
        end        
        coll(i,m)=b(1)/(f(i+1)-f(i));            
    end    
end
figure
plot(coll(:,2:8),'linewidth',2)
legend(['7654321']')