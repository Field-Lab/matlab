%% Natural Movies

clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20130302_1'
codeWord='NatMov'

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
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
    end
end

units=dir([mainpath, 'units/*.mat']);
nm=zeros(26000,size(file_list,2),length(units));
for cnt=1:length(units)
    
    load([mainpath,'units/',units(cnt).name]);
    
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            convSpikes=convolved(spikes,40,26000); 
            nm(:,i,cnt)=convSpikes(121:end-120);            
        end
    end
    if isempty(regexp(units(cnt).name,date, 'once'))
        names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
    else
        names{cnt}=units(cnt).name(1:end-4);
    end
end

save([path2save,'natMov'],'nm','names')

clear
date='20130302_1'
mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'natMov'])

namesNM=names;
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all','onOff','names','bf','wf')
for j=1:length(namesNM)
    for i=1:length(names)
        if ~isempty(regexp(namesNM{j},names{i}))
            cnt_name(j)=i;
            onOff_loc(j)=onOff(i);
            break
        end
    end
end

cc=1;
for i=1:6:size(nm,2)
    nm_av(:,cc,:)=mean(nm(:,i:i+5,:),2);
    cc=cc+1;
end


path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/NM/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end
nds='76543'
figure
set(gcf,'position',[20 50        1639         900])
for j=1:length(onOff_loc)
    cc=1;
    for i=4:3:16
        subplot('Position',[0.02 0.02+(0.2*(cc-1)) 0.68 0.16])
        plot(nm_av(:,i:i+2,j))
        if onOff_loc(j)>0
            tit='ON';
            col='r';
        elseif onOff_loc(j)<0
            tit='OFF';
            col='b';
        else
            tit='n/a';
            col='g';
        end
        title(['ND', nds(cc),'   ',int2str(cnt_name(j)),'  ',names{cnt_name(j)},'  ',tit],'Interpreter','none')
        axis tight
        subplot('Position',[0.74 0.02+(0.2*(cc-1)) 0.23 0.16])
        a=wf{cc+1}(:,cnt_name(j));
        b=bf{cc+1}(:,cnt_name(j));
        a=[a; zeros(200,1); b];
        plot(a,'color',col,'linewidth',2)
        line([0 9200],[0 0],'color','k')
        axis tight
        cc=cc+1;
    end
    saveas(gcf,[path2save,names{cnt_name(j)},'.png'])
end


