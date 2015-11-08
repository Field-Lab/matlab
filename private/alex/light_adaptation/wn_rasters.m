

datarun = load_data(fullfile(server_path(), '2015-02-24-5/data003/data003'));
datarun = load_params(datarun,'verbose',1);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);
datarun = load_ei(datarun, 'all');

datarun2 = load_data(fullfile(server_path(), '2015-02-24-5/data008/data008'));
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = set_polarities(datarun2);
datarun2 = load_neurons(datarun2);
datarun2 = load_ei(datarun2, 'all');


datarun3 = load_data(fullfile(server_path(), '2015-02-24-5/data014/data014'));
datarun3 = load_params(datarun3,'verbose',1);
datarun3 = set_polarities(datarun3);
datarun3 = load_neurons(datarun3);
datarun3 = load_ei(datarun3, 'all');

datarun4 = load_data(fullfile(server_path(), '2015-02-24-5/data017/data017'));
datarun4 = load_params(datarun4,'verbose',1);
datarun4 = set_polarities(datarun4);
datarun4 = load_neurons(datarun4);
datarun4 = load_ei(datarun4, 'all');

a = map_ei(datarun, datarun2);
b = map_ei(datarun, datarun3);
c = map_ei(datarun, datarun4);


masterTrigs=datarun.triggers;
figure
plot(diff(masterTrigs))

myTrigs=find(diff(masterTrigs)>0.84&diff(masterTrigs)<1);
myTrigs=myTrigs(myTrigs<740);
myTrigs=[0; myTrigs];

slaveTrigs=datarun2.triggers;
figure
plot(diff(slaveTrigs))

mySTrigs=find(diff(slaveTrigs)>0.84&diff(slaveTrigs)<1);
mySTrigs=mySTrigs(mySTrigs<740);
mySTrigs=[0; mySTrigs];

file_path='/Volumes/Analysis/2015-02-24-5/wn_rasters/data003_008_014_017/';
if ~exist(file_path,'dir')
    mkdir(file_path);
end

for i=1:length(a)
    
    if ~isempty(a{i})&&~isempty(b{i})&&~isempty(c{i})

        masterCell=i;
        masterSpikes=datarun.spikes{masterCell};
        
        slaveCell=find(datarun2.cell_ids==a{i});
        slaveSpikes=datarun2.spikes{slaveCell};
        
        slaveCell3=find(datarun3.cell_ids==b{i});
        slaveSpikes3=datarun3.spikes{slaveCell3};
        
        slaveCell4=find(datarun4.cell_ids==c{i});
        slaveSpikes4=datarun4.spikes{slaveCell4};
        
        
        myMasterSpikes=cell(19,1);
        mySlaveSpikes=cell(19,1);
        mySlaveSpikes3=cell(19,1);
        mySlaveSpikes4=cell(19,1);
        
        for j=1:19
            tmp=masterSpikes(masterSpikes>masterTrigs(myTrigs(j)+1) & masterSpikes<masterTrigs(myTrigs(j+1)))...
                - masterTrigs(myTrigs(j)+1);
            myMasterSpikes{j}=tmp;
            
            tmp=slaveSpikes(slaveSpikes>slaveTrigs(mySTrigs(j)+1) & slaveSpikes<slaveTrigs(mySTrigs(j+1)))...
                - slaveTrigs(mySTrigs(j)+1);
            mySlaveSpikes{j}=tmp;
            
            tmp=slaveSpikes3(slaveSpikes3>slaveTrigs(mySTrigs(j)+1) & slaveSpikes3<slaveTrigs(mySTrigs(j+1)))...
                - slaveTrigs(mySTrigs(j)+1);
            mySlaveSpikes3{j}=tmp;
            
            tmp=slaveSpikes4(slaveSpikes4>slaveTrigs(mySTrigs(j)+1) & slaveSpikes4<slaveTrigs(mySTrigs(j+1)))...
                - slaveTrigs(mySTrigs(j)+1);
            mySlaveSpikes4{j}=tmp;
        end

        t=[];
        cnt=0;fr=0;fr1=0;fr3=0;fr4=0;
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
        for j=1:19
            t=[t mySlaveSpikes3{j}'*1000+31000*cnt];
            cnt=cnt+1;
            tmp=convolved(mySlaveSpikes3{j}'*1000,40,31000);
            fr3=fr3+tmp(((size(tmp,2)-31000)/2+1):end-((size(tmp,2)-31000)/2));
        end
        for j=1:19
            t=[t mySlaveSpikes4{j}'*1000+31000*cnt];
            cnt=cnt+1;
            tmp=convolved(mySlaveSpikes4{j}'*1000,40,31000);
            fr4=fr4+tmp(((size(tmp,2)-31000)/2+1):end-((size(tmp,2)-31000)/2));
        end

        
        figure
        set(gcf,'position',[82         248        1158         850])
        h=subplot(2,1,1);
        rasterplot(t,cnt,31000,h)
        line([0,31000], [28.25,28.25],'color','r','linewidth',1.2)
        line([0,31000], [56.75,56.75],'color','r','linewidth',1.2)
        line([0,31000], [85.25,85.25],'color','r','linewidth',1.2)
        set(gca,'ytick',14:28:98,'yticklabel',{'NDF3.0','NDF2.0'})
        set(gca,'xtick',0:5000:30000,'xticklabel',{'0','5','10','15','20','25','30'})
        ylabel('light level')
        xlabel('NSEM,s')
        title(['2015-02-24-5, data003 cell ',int2str(datarun.cell_ids(i)),', data008 cell ',int2str(a{i})])
        
        subplot(2,1,2)
        hold on
        plot(fr/20,'r','linewidth',1.2)
        plot(fr1/20,'b','linewidth',1.2)
        plot(fr3/20,'m','linewidth',1.2)
        plot(fr4/20,'c','linewidth',1.2)
        axis([0 31000 0 Inf])
        set(gca,'xtick',0:5000:30000,'xticklabel',{'0','5','10','15','20','25','30'})
        legend('NDF 3.0', 'NDF 2.0');

        
        saveas(gcf,[file_path, 'data003_cell',int2str(datarun.cell_ids(i)),'.jpg'])
        close(gcf)
    end
end

