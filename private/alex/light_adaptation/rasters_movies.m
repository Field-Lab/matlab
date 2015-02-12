
datarun = load_data(fullfile(server_path(), '2015-01-29-0/data002/data002'));
datarun = load_params(datarun,'verbose',1);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);
datarun = load_ei(datarun, 'all');

datarun2 = load_data(fullfile(server_path(), '2015-01-29-0/data004/data004'));
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = set_polarities(datarun2);
datarun2 = load_neurons(datarun2);
datarun2 = load_ei(datarun2, 'all');

a = map_ei(datarun, datarun2);

masterTrigs=datarun.triggers;
figure
plot(diff(masterTrigs))

myTrigs=find(diff(masterTrigs)>0.9&diff(masterTrigs)<2);
myTrigs=myTrigs(myTrigs<740);
myTrigs=[0; myTrigs];

slaveTrigs=datarun2.triggers;

mySTrigs=find(diff(slaveTrigs)>0.9&diff(slaveTrigs)<2);
mySTrigs=mySTrigs(mySTrigs<740);
mySTrigs=[0; mySTrigs];

for i=1:length(a)
    
    if ~isempty(a{i})
        masterCell=i;
        slaveCell=find(datarun2.cell_ids==a{i});
        masterSpikes=datarun.spikes{masterCell};
        slaveSpikes=datarun2.spikes{slaveCell};
        myMasterSpikes=cell(19,1);
        mySlaveSpikes=cell(19,1);
        for j=1:19
            tmp=masterSpikes(masterSpikes>masterTrigs(myTrigs(j)+1) & masterSpikes<masterTrigs(myTrigs(j+1)))...
                - masterTrigs(myTrigs(j)+1);
            myMasterSpikes{j}=tmp;
            
            tmp=slaveSpikes(slaveSpikes>slaveTrigs(mySTrigs(j)+1) & slaveSpikes<slaveTrigs(mySTrigs(j+1)))...
                - slaveTrigs(mySTrigs(j)+1);
            mySlaveSpikes{j}=tmp;
        end

        t=[];
        cnt=0;
        fr=0;fr1=0;
        for j=1:19
            t=[t myMasterSpikes{j}'*1000+31000*cnt];
            cnt=cnt+1;
            tmp=convolved(myMasterSpikes{j}'*1000,40,31000);
            fr=fr+tmp(((size(tmp,2)-31000)/2+1):end-((size(tmp,2)-31000)/2));
        end
        for j=1:19
            t=[t mySlaveSpikes{j}'*1000+31000*cnt];
            cnt=cnt+1;
            tmp=convolved(mySlaveSpikes{j}'*1000,40,31000);
            fr1=fr1+tmp(((size(tmp,2)-31000)/2+1):end-((size(tmp,2)-31000)/2));
        end
%         rasterplot(t,cnt,31000)
%         line([0,31000], [28.25,28.25],'color','r','linewidth',2)
%         line([0,31000], [28.25,28.25],'color','r','linewidth',2)
%         set(gca,'ytick',14:28:43,'yticklabel',{'NDF3.0','NDF2.0'},'fontsize',16,'fontweight','bold')
%         set(gca,'xtick',0:5000:30000,'xticklabel',{'0','5','10','15','20','25','30'},'fontsize',16,'fontweight','bold')
%         ylabel('light level')
%         xlabel('NSEM,s')
%         title(['2015-01-29-0, data002 cell ',int2str(datarun.cell_ids(i)),', data004 cell ',int2str(a{i})])
%         saveas(gcf,['/Users/alexth/Desktop/Light_adaptation/movie_rasters/2015-01-29-0/data002_cell',int2str(datarun.cell_ids(i)),'.jpg'])
%         close(gcf)
%         
%         figure
%         set(gcf,'position',[31         717        1810         381])
%         hold on
%         plot(fr/20,'r','linewidth',1.2)
%         plot(fr1/20,'b','linewidth',1.2)
%         axis([0 30000 0 Inf])
%         legend('NDF 3.0', 'NDF 2.0');
%         title(['2015-01-29-0, data002 cell ',int2str(datarun2.cell_ids(i)),', data004 cell ',int2str(a{i})])
%         saveas(gcf,['/Users/alexth/Desktop/Light_adaptation/movie_fr/2015-01-29-0/data002_cell',int2str(datarun.cell_ids(i)),'.jpg'])
%         close(gcf)

        figure
        set(gcf,'position',[82         248        1158         850])
        h=subplot(2,1,1);
        rasterplot(t,cnt,31000,h)
        line([0,31000], [28.25,28.25],'color','r','linewidth',1.2)
        set(gca,'ytick',14:28:43,'yticklabel',{'NDF3.0','NDF2.0'})
        set(gca,'xtick',0:5000:30000,'xticklabel',{'0','5','10','15','20','25','30'})
        ylabel('light level')
        xlabel('NSEM,s')
        title(['2015-01-29-0, data002 cell ',int2str(datarun.cell_ids(i)),', data004 cell ',int2str(a{i})])
        
        subplot(2,1,2)
        hold on
        plot(fr/20,'r','linewidth',1.2)
        plot(fr1/20,'b','linewidth',1.2)
        axis([0 30000 0 Inf])
        set(gca,'xtick',0:5000:30000,'xticklabel',{'0','5','10','15','20','25','30'})
        legend('NDF 3.0', 'NDF 2.0');

        
        saveas(gcf,['/Users/alexth/Desktop/Light_adaptation/movie_rasters_fr_comb/2015-01-29-0/data002_cell',int2str(datarun.cell_ids(i)),'.jpg'])
        close(gcf)
    end
end


datarun3 = load_data(fullfile(server_path(), '2015-01-29-0/data001/data001'));
datarun3 = load_sta(datarun3,'load_sta',[],'keep_java_sta',true);
datarun3 = load_params(datarun3,'verbose',1);
datarun3 = set_polarities(datarun3);
datarun3 = load_neurons(datarun3);
datarun3 = load_ei(datarun3, 'all');

b = map_ei(datarun, datarun3);

rgcs=[16 3817 4804 5329 7128];
cellT=[];
for i=1:length(rgcs)
    tmp=find(datarun.cell_ids==rgcs(i));
    newCell=b{tmp}
    if ~isempty(newCell)
        newtmp=find(datarun3.cell_ids==newCell);
        for j=1:10
            d=find(datarun3.cell_types{j}.cell_ids==newCell);
            if ~isempty(d)
                cellT=[cellT j];
                break;
            end
        end
    else
        cellT=[cellT 0];
    end
end


datarun4 = load_data(fullfile(server_path(), '2015-01-29-0/data003/data003'));
datarun4 = load_sta(datarun4,'load_sta',[],'keep_java_sta',true);
datarun4 = load_params(datarun4,'verbose',1);
datarun4 = set_polarities(datarun4);
datarun4 = load_neurons(datarun4);
datarun4 = load_ei(datarun4, 'all');

c = map_ei(datarun, datarun4);


rgcs=[16 3817 4804 5329 7128];
cellT=[];
for i=1:length(rgcs)
    tmp=find(datarun.cell_ids==rgcs(i));
    newCell=c{tmp}
    if ~isempty(newCell)
        newtmp=find(datarun4.cell_ids==newCell);
        for j=1:10
            d=find(datarun4.cell_types{j}.cell_ids==newCell);
            if ~isempty(d)
                cellT=[cellT j];
                break;
            end
        end
    else
        cellT=[cellT 0];
    end
end






plot_rf_summaries(datarun3, [286], 'clear', false,  'plot_fits', true, 'fit_color', 'r')
hold on
plot_rf_summaries(datarun4, [286], 'clear', false,  'plot_fits', true, 'fit_color', 'b')
legend('NDF 3', 'NDF 2')

figure
plot_rf_summaries(datarun3, {1}, 'clear', false,  'plot_fits', true, 'fit_color', 'r')
hold on
plot_rf_summaries(datarun4, {1}, 'clear', false,  'plot_fits', true, 'fit_color', 'b')
legend('NDF 3', 'NDF 2')


%****************


datarun = load_data(fullfile(server_path(), '2015-01-29-2/data001/data001'));
datarun = load_params(datarun,'verbose',1);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);
datarun = load_ei(datarun, 'all');

datarun2 = load_data(fullfile(server_path(), '2015-01-29-2/data003/data003'));
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = set_polarities(datarun2);
datarun2 = load_neurons(datarun2);
datarun2 = load_ei(datarun2, 'all');

a = map_ei(datarun, datarun2);

masterTrigs=datarun.triggers;
figure
plot(diff(masterTrigs))

myTrigs=find(diff(masterTrigs)>0.9&diff(masterTrigs)<2);
myTrigs=myTrigs(myTrigs<740);
myTrigs=[0; myTrigs];

slaveTrigs=datarun2.triggers;

mySTrigs=find(diff(slaveTrigs)>0.9&diff(slaveTrigs)<2);
mySTrigs=mySTrigs(mySTrigs<740);
mySTrigs=[0; mySTrigs];

for i=1:length(a)
    
    if ~isempty(a{i})
        masterCell=i;
        slaveCell=find(datarun2.cell_ids==a{i});
        masterSpikes=datarun.spikes{masterCell};
        slaveSpikes=datarun2.spikes{slaveCell};
        myMasterSpikes=cell(19,1);
        mySlaveSpikes=cell(19,1);
        for j=1:19
            tmp=masterSpikes(masterSpikes>masterTrigs(myTrigs(j)+1) & masterSpikes<masterTrigs(myTrigs(j+1)))...
                - masterTrigs(myTrigs(j)+1);
            myMasterSpikes{j}=tmp;
            
            tmp=slaveSpikes(slaveSpikes>slaveTrigs(mySTrigs(j)+1) & slaveSpikes<slaveTrigs(mySTrigs(j+1)))...
                - slaveTrigs(mySTrigs(j)+1);
            mySlaveSpikes{j}=tmp;
        end

        t=[];
        cnt=0;fr=0;fr1=0;
        for j=1:19
            t=[t myMasterSpikes{j}'*1000+31000*cnt];
            cnt=cnt+1;
            tmp=convolved(myMasterSpikes{j}'*1000,40,31000);
            fr=fr+tmp(((size(tmp,2)-31000)/2+1):end-((size(tmp,2)-31000)/2));
        end
        for j=1:19
            t=[t mySlaveSpikes{j}'*1000+31000*cnt];
            cnt=cnt+1;
            tmp=convolved(mySlaveSpikes{j}'*1000,40,31000);
            fr1=fr1+tmp(((size(tmp,2)-31000)/2+1):end-((size(tmp,2)-31000)/2));
        end
%         rasterplot(t,cnt,31000)
%         line([0,31000], [28.25,28.25],'color','r','linewidth',2)
%         line([0,31000], [28.25,28.25],'color','r','linewidth',2)
%         set(gca,'ytick',14:28:43,'yticklabel',{'NDF3.0','NDF2.0'},'fontsize',16,'fontweight','bold')
%         set(gca,'xtick',0:5000:30000,'xticklabel',{'0','5','10','15','20','25','30'},'fontsize',16,'fontweight','bold')
%         ylabel('light level')
%         xlabel('NSEM,s')
%         title(['2015-01-29-2, data001 cell ',int2str(datarun.cell_ids(i)),', data003 cell ',int2str(a{i})])
%         saveas(gcf,['/Users/alexth/Desktop/Light_adaptation/movie_rasters/2015-01-29-2/data001_cell',int2str(datarun.cell_ids(i)),'.jpg'])
%         close(gcf)
        
%         figure
%         set(gcf,'position',[31         717        1810         381])
%         hold on
%         plot(fr/20,'r','linewidth',1.2)
%         plot(fr1/20,'b','linewidth',1.2)
%         axis([0 30000 0 Inf])
%         legend('NDF 3.0', 'NDF 2.0');
%         title(['2015-01-29-2, data001 cell ',int2str(datarun.cell_ids(i)),', data003 cell ',int2str(a{i})])
%         saveas(gcf,['/Users/alexth/Desktop/Light_adaptation/movie_fr/2015-01-29-2/data001_cell',int2str(datarun.cell_ids(i)),'.jpg'])
%         close(gcf)
        
        figure
        set(gcf,'position',[82         248        1158         850])
        h=subplot(2,1,1);
        rasterplot(t,cnt,31000,h)
        line([0,31000], [28.25,28.25],'color','r','linewidth',1.2)
        set(gca,'ytick',14:28:43,'yticklabel',{'NDF3.0','NDF2.0'})
        set(gca,'xtick',0:5000:30000,'xticklabel',{'0','5','10','15','20','25','30'})
        ylabel('light level')
        xlabel('NSEM,s')
        title(['2015-01-29-2, data001 cell ',int2str(datarun.cell_ids(i)),', data003 cell ',int2str(a{i})])
        
        subplot(2,1,2)
        hold on
        plot(fr/20,'r','linewidth',1.2)
        plot(fr1/20,'b','linewidth',1.2)
        axis([0 30000 0 Inf])
        set(gca,'xtick',0:5000:30000,'xticklabel',{'0','5','10','15','20','25','30'})
        legend('NDF 3.0', 'NDF 2.0');

        
        saveas(gcf,['/Users/alexth/Desktop/Light_adaptation/movie_rasters_fr_comb/2015-01-29-2/data001_cell',int2str(datarun.cell_ids(i)),'.jpg'])
        close(gcf)
    end
end

