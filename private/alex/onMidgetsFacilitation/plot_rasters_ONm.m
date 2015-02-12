date='2011-12-13-2';

% STA
run='data004-0';
path2load = fullfile(server_path(), [date, '/streamed/',run,'/',run]);
datarun2 = load_data(path2load);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2,'load_sta',[],'keep_java_sta',true);
datarun2 = set_polarities(datarun2);
datarun2 = load_ei(datarun2, 'all');


% UDCR
run='data007';
path2load = fullfile(server_path(), [date, '/',run,'/',run]);
datarun = load_data(path2load);
datarun = load_params(datarun,'verbose',1);
datarun = set_polarities(datarun);
datarun = load_ei(datarun, 'all');
datarun = load_neurons(datarun);

stimulus=read_stim_lisp_output_ath('2011-12-13-2','s07');
parsed=parse_stim_rgbs_ath(stimulus);
map=load(parsed.mappath);
figure
imagesc(map)
% find stimuli patterns
myStim=zeros(4,length(stimulus.rgbs));
for i=1:length(stimulus.rgbs) 
    myStim(1:4,i)=stimulus.rgbs{i}(:,1);
end

[B, ic, ia] = unique(myStim', 'rows'); % 73 patterns, 40 repetitions

a = map_ei(datarun, datarun2);



cellID=6601;
myCellSTA=find(datarun2.cell_ids==cellID);
for i=1:length(a)
    if a{i}==cellID
        myCell=i;        
        break
    end
end


for myCell=271:280
    % calculate response
    spikes=datarun.spikes{myCell}*1000;
    tr=datarun.triggers(1:2:end)*1000;
    
    mySpikes=[];
    for i=1:length(ic)
        patts=find(ia==i);
        cnt=1;
        for j=1:length(patts)
            mySpikes{cnt,i}=spikes(spikes>tr(patts(j))&spikes<tr(patts(j))+750)-tr(patts(j));
            cnt=cnt+1;
        end
    end
    
    figure
    for i=1:4
        for j=1:4
            if i==j
                patN=find(B(:,i)==0.48&sum(B,2)==0.48);
                patN2=find(B(:,i)==0.24&sum(B,2)==0.24);
                patN3=find(B(:,i)==-0.48&sum(B,2)==-0.48);
            else
                patN=find(B(:,i)==0.48&B(:,j)==-0.48);
                patN2=find(B(:,i)==0.24&B(:,j)==-0.48);
                patN3=find(B(:,i)==0.12&B(:,j)==-0.48);
            end
            
            
            t=[];
            cnt=0;
            fr=0;
            fr2=0;fr3=0;
            for k=1:40
                t=[t mySpikes{k,patN}'+800*cnt];
                
                tmp=convolved(mySpikes{k,patN},10,800);
                fr=fr+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                tmp=convolved(mySpikes{k,patN2},10,800);
                fr2=fr2+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                tmp=convolved(mySpikes{k,patN3},10,800);
                fr3=fr3+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                
                cnt=cnt+1;
            end
            sbind=sub2ind([4,4],j,i);
            h=subplot(4,4,sbind);
            hold on
            plot(fr/20,'r','linewidth',2)
            plot(fr2/20,'b','linewidth',2)
            plot(fr3/20,'g','linewidth',2)
            if sbind==1
                legend('1','0.5','-1')
            elseif sbind==2
                legend('1,-1','0.5,-1','0.25,-1')
            end
            rasterplot(t,cnt,800,h)
            
        end
    end
end


respcells=[9 12 13 16 17 19 22 23 24 29 30 31 43 47 48 51 55 58 59 61  62 63 ...
    66 67 70 72 73 74 78 81 82 83 84 86 87 89 91 93 95 96 98 99 100 103 108 109 110 ...
    115 117 119 121 122 128 130 137 138 147 169 171 177 181 187 194 200 205 207 208 ...
    211 219 221 223 225 227 228 230 232 233 235 236 237 241 243 244 246 247 249 ...
    251 255 259 260 261 274];


for myCell=respcells
    % calculate response
    spikes=datarun.spikes{myCell}*1000;
    tr=datarun.triggers(1:2:end)*1000;
    
    cellMap(myCell) % cell ID in STA datarun
    
    tmp=find(cellTypes(1,:)==cellMap(myCell));
    
    if ~isempty(tmp)
        myCellType=['STA ID ',int2str(tmp),', ',datarun2.cell_types{cellTypes(2,tmp)}.name];
    else
        myCellType='no slave';
    end
    tit=[int2str(myCell),', ID ',int2str(datarun.cell_ids(myCell)), ', ',myCellType];
    
    mySpikes=[];
    for i=1:length(ic)
        patts=find(ia==i);
        cnt=1;
        for j=1:length(patts)
            mySpikes{cnt,i}=spikes(spikes>tr(patts(j))&spikes<tr(patts(j))+750)-tr(patts(j));
            cnt=cnt+1;
        end
    end
    
    figure
    set(gcf,'position',[1           1        1920        1105])
    for i=1:4
        for j=1:4
            if i==j
                patN=find(B(:,i)==0.48&sum(B,2)==0.48);
                patN2=find(B(:,i)==0.24&sum(B,2)==0.24);
                patN3=find(B(:,i)==-0.48&sum(B,2)==-0.48);
            else
                patN=find(B(:,i)==0.48&B(:,j)==-0.48);
                patN2=find(B(:,i)==0.24&B(:,j)==-0.48);
                patN3=find(B(:,i)==0.12&B(:,j)==-0.48);
            end
            
            
            t=[];
            cnt=0;
            fr=0;
            fr2=0;fr3=0;
            for k=1:40
                t=[t mySpikes{k,patN}'+800*cnt];
                
                tmp=convolved(mySpikes{k,patN},10,800);
                fr=fr+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                tmp=convolved(mySpikes{k,patN2},10,800);
                fr2=fr2+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                tmp=convolved(mySpikes{k,patN3},10,800);
                fr3=fr3+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                
                cnt=cnt+1;
            end
            sbind=sub2ind([4,4],j,i);
            h=subplot(4,4,sbind);
            hold on
            plot(fr/20,'r','linewidth',2)
            plot(fr2/20,'b','linewidth',2)
            plot(fr3/20,'g','linewidth',2)
            if sbind==1
                hl=legend('1','0.5','-1');
                set(hl, 'Box', 'off')
                set(hl, 'Color', 'none')
                set(hl,'location','best')
                title(tit)
            elseif sbind==2
                hl=legend('1,-1','0.5,-1','0.25,-1');
                set(hl, 'Box', 'off')
                set(hl, 'Color', 'none')
                set(hl,'location','best')
            end
            rasterplot(t,cnt,800,h)
            
        end
    end
    
    saveas(gcf,['/Users/alexth/Desktop/ONmidgets_facilitation/',date,'/cell_', ...
        int2str(myCell),'_ID_',int2str(datarun.cell_ids(myCell)),'.jpg'])
    close all
    
end

cellTypes=zeros(2,length(datarun2.cell_ids));
cellTypes(1,:)=datarun2.cell_ids;
for j=1:length(datarun2.cell_types)
    if ~isempty(datarun2.cell_types{j}.cell_ids)
        [~,ia, ib]=intersect(cellTypes(1,:),datarun2.cell_types{j}.cell_ids);
        cellTypes(2,ia)=j;
    end
end

cellMap=zeros(size(datarun.cell_ids));cnt=1;
for myCell=1:length(datarun.cell_ids)
    if isempty(a{myCell})
        cellMap(cnt)=0;
    else
        cellMap(cnt)=a{myCell};
        find(datarun2.cell_ids==cellMap(cnt));
    end
    cnt=cnt+1;
end
    
[~,ia,ib]=intersect(cellTypes(1,:),cellMap);
cellTypes(:,ia);



%%
date='2011-12-13-2';

% STA
run='data004_data007_data008-norefit';
tmp='/Volumes/Analysis/2011-12-13-2/data004_data007_data008-norefit/data004-from-data004_data007_data008/data004-from-data004_data007_data008'
datarun2 = load_data(tmp);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2,'load_sta',[],'keep_java_sta',true);
datarun2 = set_polarities(datarun2);
datarun2 = load_ei(datarun2, 'all');


% UDCR
run='data004_data007_data008-norefit';
tmp='/Volumes/Analysis/2011-12-13-2/data004_data007_data008-norefit/data007-from-data004_data007_data008/data007-from-data004_data007_data008'
datarun = load_data(tmp);
datarun = load_params(datarun,'verbose',1);
datarun = set_polarities(datarun);
datarun = load_ei(datarun, 'all');
datarun = load_neurons(datarun);

stimulus=read_stim_lisp_output_ath('2011-12-13-2','s07');
parsed=parse_stim_rgbs_ath(stimulus);
map=load(parsed.mappath);
figure
imagesc(map)
% find stimuli patterns
myStim=zeros(4,length(stimulus.rgbs));
for i=1:length(stimulus.rgbs) 
    myStim(1:4,i)=stimulus.rgbs{i}(:,1);
end

[B, ic, ipat] = unique(myStim', 'rows'); % 73 patterns, 40 repetitions


cellTypes=zeros(1,length(datarun2.cell_ids));
for j=1:length(datarun2.cell_types)
    if ~isempty(datarun2.cell_types{j}.cell_ids)
        [~,ia, ib]=intersect(datarun2.cell_ids,datarun2.cell_types{j}.cell_ids);
        cellTypes(ia)=j;
    end
end

for myCell=1:length(datarun.cell_ids)
    % calculate response
    spikes=datarun.spikes{myCell}*1000;
    tr=datarun.triggers(1:2:end)*1000;
            
    myCellType=datarun2.cell_types{cellTypes(myCell)}.name;

    tit=[int2str(myCell),', ID ',int2str(datarun.cell_ids(myCell)), ', ',myCellType];
    
    mySpikes=[];
    for i=1:length(ic)
        patts=find(ipat==i);
        cnt=1;
        for j=1:length(patts)
            mySpikes{cnt,i}=spikes(spikes>tr(patts(j))&spikes<tr(patts(j))+750)-tr(patts(j));
            cnt=cnt+1;
        end
    end
    
    figure
    set(gcf,'position',[1           1        1920        1105])
    for i=1:4
        for j=1:4
            if i==j
                patN=find(B(:,i)==0.48&sum(B,2)==0.48);
                patN2=find(B(:,i)==0.24&sum(B,2)==0.24);
                patN3=find(B(:,i)==-0.48&sum(B,2)==-0.48);
            else
                patN=find(B(:,i)==0.48&B(:,j)==-0.48);
                patN2=find(B(:,i)==0.24&B(:,j)==-0.48);
                patN3=find(B(:,i)==0.12&B(:,j)==-0.48);
            end
            
            
            t=[];
            cnt=0;
            fr=0;
            fr2=0;fr3=0;
            for k=1:40
                t=[t mySpikes{k,patN}'+800*cnt];
                
                tmp=convolved(mySpikes{k,patN},10,800);
                fr=fr+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                tmp=convolved(mySpikes{k,patN2},10,800);
                fr2=fr2+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                tmp=convolved(mySpikes{k,patN3},10,800);
                fr3=fr3+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                
                cnt=cnt+1;
            end
            sbind=sub2ind([4,4],j,i);
            h=subplot(4,4,sbind);
            hold on
            plot(fr/20,'r','linewidth',2)
            plot(fr2/20,'b','linewidth',2)
            plot(fr3/20,'g','linewidth',2)
            if sbind==1
                hl=legend('1','0.5','-1');
                set(hl, 'Box', 'off')
                set(hl, 'Color', 'none')
                set(hl,'location','best')
                title(tit)
            elseif sbind==2
                hl=legend('1,-1','0.5,-1','0.25,-1');
                set(hl, 'Box', 'off')
                set(hl, 'Color', 'none')
                set(hl,'location','best')
            end
            rasterplot(t,cnt,800,h)
            
        end
    end
    
    saveas(gcf,['/Users/alexth/Desktop/ONmidgets_facilitation/',date,'/cell_', ...
        int2str(myCell),'_ID_',int2str(datarun.cell_ids(myCell)),'.jpg'])
    close all
    
end


%% 2011-12-13-2 OFF midgets
date='2011-12-13-2';

% STA
run='data008-from-data008_data012_data013_data014';
tmp='/Volumes/Analysis/2011-12-13-2/data008_data012_data013_data014-norefit/data008-from-data008_data012_data013_data014/data008-from-data008_data012_data013_data014'
datarun2 = load_data(tmp);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2,'load_sta',[],'keep_java_sta',true);
datarun2 = set_polarities(datarun2);
datarun2 = load_ei(datarun2, 'all');


% UDCR
run='data013-from-data008_data012_data013_data014';
tmp='/Volumes/Analysis/2011-12-13-2/data008_data012_data013_data014-norefit/data013-from-data008_data012_data013_data014/data013-from-data008_data012_data013_data014'
datarun = load_data(tmp);
datarun = load_params(datarun,'verbose',1);
datarun = set_polarities(datarun);
datarun = load_ei(datarun, 'all');
datarun = load_neurons(datarun);

stimulus=read_stim_lisp_output_ath('2011-12-13-2','s13');
parsed=parse_stim_rgbs_ath(stimulus);
map=load(parsed.mappath);
figure
imagesc(map)
% find stimuli patterns
myStim=zeros(4,length(stimulus.rgbs));
for i=1:length(stimulus.rgbs) 
    myStim(1:4,i)=stimulus.rgbs{i}(:,1);
end

[B, ic, ipat] = unique(myStim', 'rows'); % 17 patterns, 40 repetitions


cellTypes=zeros(1,length(datarun2.cell_ids));
for j=1:length(datarun2.cell_types)
    if ~isempty(datarun2.cell_types{j}.cell_ids)
        [~,ia, ib]=intersect(datarun2.cell_ids,datarun2.cell_types{j}.cell_ids);
        cellTypes(ia)=j;
    end
end

subplots=1;
figure
set(gcf,'position',[1           1        1920        1105])
cntcells=1;
for myCell=1:length(datarun.cell_ids)
    % calculate response
    spikes=datarun.spikes{myCell}*1000;
    tr=datarun.triggers(1:2:end)*1000;
            
    myCellType=datarun2.cell_types{cellTypes(myCell)}.name;

    tit=[int2str(myCell),', ID ',int2str(datarun.cell_ids(myCell)), ', ',myCellType];
    
    mySpikes=[];
    for i=1:length(ic)
        patts=find(ipat==i);
        cnt=1;
        for j=1:length(patts)-2
            mySpikes{cnt,i}=spikes(spikes>tr(patts(j))&spikes<tr(patts(j))+750)-tr(patts(j));
            cnt=cnt+1;
        end
    end

    for i=1:4
        patN=find(B(:,i)==-0.48);
        patN2=find(B(:,i)==B(2,1));
        patN3=find(B(:,i)==-0.24);
        patN4=find(B(:,i)==-0.12);
    end
    
    t=[];
    cnt=0;
    fr=0;
    fr2=0;fr3=0;fr4=0;
    for k=1:28
        t=[t mySpikes{k,patN}'+800*cnt];
        
        tmp=convolved(mySpikes{k,patN},10,800);
        fr=fr+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
        tmp=convolved(mySpikes{k,patN2},10,800);
        fr2=fr2+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
        tmp=convolved(mySpikes{k,patN3},10,800);
        fr3=fr3+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
        tmp=convolved(mySpikes{k,patN4},10,800);
        fr4=fr4+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
        
        cnt=cnt+1;
    end
    if subplots<17
        h=subplot(4,4,subplots);
        subplots=subplots+1;
    else
        close all
        figure
        set(gcf,'position',[1           1        1920        1105])
        subplots=1;
        h=subplot(4,4,subplots);  
        subplots=subplots+1;
    end

    hold on
    plot(fr/20,'r','linewidth',2)
    plot(fr2/20,'b','linewidth',2)
    plot(fr3/20,'g','linewidth',2)
    plot(fr4/20,'m','linewidth',2)
    title(tit)
    rasterplot(t,cnt,800,h)
    
    if subplots==17
        saveas(gcf,['/Users/alexth/Desktop/ONmidgets_facilitation/2011-12-13-2_offs_data013/cells_', ...
            int2str(cntcells),'.jpg'])
        cntcells=cntcells+16;
    end
    
  
end




%% 2012-09-13-2
date=' 2012-09-13-2';


% STA
run='data004_data007_data008-norefit';
tmp='/Volumes/Analysis/2011-12-13-2/data004_data007_data008-norefit/data004-from-data004_data007_data008/data004-from-data004_data007_data008'
datarun2 = load_data(tmp);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2,'load_sta',[],'keep_java_sta',true);
datarun2 = set_polarities(datarun2);
datarun2 = load_ei(datarun2, 'all');


% UDCR
run='data004_data007_data008-norefit';
tmp='/Volumes/Analysis/2011-12-13-2/data004_data007_data008-norefit/data007-from-data004_data007_data008/data007-from-data004_data007_data008'
datarun = load_data(tmp);
datarun = load_params(datarun,'verbose',1);
datarun = set_polarities(datarun);
datarun = load_ei(datarun, 'all');
datarun = load_neurons(datarun);

stimulus=read_stim_lisp_output_ath('2011-12-13-2','s07');
parsed=parse_stim_rgbs_ath(stimulus);
map=load(parsed.mappath);
figure
imagesc(map)
% find stimuli patterns
myStim=zeros(4,length(stimulus.rgbs));
for i=1:length(stimulus.rgbs) 
    myStim(1:4,i)=stimulus.rgbs{i}(:,1);
end

[B, ic, ipat] = unique(myStim', 'rows'); % 73 patterns, 40 repetitions


cellTypes=zeros(1,length(datarun2.cell_ids));
for j=1:length(datarun2.cell_types)
    if ~isempty(datarun2.cell_types{j}.cell_ids)
        [~,ia, ib]=intersect(datarun2.cell_ids,datarun2.cell_types{j}.cell_ids);
        cellTypes(ia)=j;
    end
end

for myCell=1:length(datarun.cell_ids)
    % calculate response
    spikes=datarun.spikes{myCell}*1000;
    tr=datarun.triggers(1:2:end)*1000;
            
    myCellType=datarun2.cell_types{cellTypes(myCell)}.name;

    tit=[int2str(myCell),', ID ',int2str(datarun.cell_ids(myCell)), ', ',myCellType];
    
    mySpikes=[];
    for i=1:length(ic)
        patts=find(ipat==i);
        cnt=1;
        for j=1:length(patts)
            mySpikes{cnt,i}=spikes(spikes>tr(patts(j))&spikes<tr(patts(j))+750)-tr(patts(j));
            cnt=cnt+1;
        end
    end
    
    figure
    set(gcf,'position',[1           1        1920        1105])
    for i=1:4
        for j=1:4
            if i==j
                patN=find(B(:,i)==0.48&sum(B,2)==0.48);
                patN2=find(B(:,i)==0.24&sum(B,2)==0.24);
                patN3=find(B(:,i)==-0.48&sum(B,2)==-0.48);
            else
                patN=find(B(:,i)==0.48&B(:,j)==-0.48);
                patN2=find(B(:,i)==0.24&B(:,j)==-0.48);
                patN3=find(B(:,i)==0.12&B(:,j)==-0.48);
            end
            
            
            t=[];
            cnt=0;
            fr=0;
            fr2=0;fr3=0;
            for k=1:40
                t=[t mySpikes{k,patN}'+800*cnt];
                
                tmp=convolved(mySpikes{k,patN},10,800);
                fr=fr+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                tmp=convolved(mySpikes{k,patN2},10,800);
                fr2=fr2+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                tmp=convolved(mySpikes{k,patN3},10,800);
                fr3=fr3+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                
                cnt=cnt+1;
            end
            sbind=sub2ind([4,4],j,i);
            h=subplot(4,4,sbind);
            hold on
            plot(fr/20,'r','linewidth',2)
            plot(fr2/20,'b','linewidth',2)
            plot(fr3/20,'g','linewidth',2)
            if sbind==1
                hl=legend('1','0.5','-1');
                set(hl, 'Box', 'off')
                set(hl, 'Color', 'none')
                set(hl,'location','best')
                title(tit)
            elseif sbind==2
                hl=legend('1,-1','0.5,-1','0.25,-1');
                set(hl, 'Box', 'off')
                set(hl, 'Color', 'none')
                set(hl,'location','best')
            end
            rasterplot(t,cnt,800,h)
            
        end
    end
    
    saveas(gcf,['/Users/alexth/Desktop/ONmidgets_facilitation/',date,'/cell_', ...
        int2str(myCell),'_ID_',int2str(datarun.cell_ids(myCell)),'.jpg'])
    close all
    
end


%% 2012-09-24-5
date=' 2012-09-24-5';


% STA
tmp='/Volumes/Analysis/2012-09-24-5/d01-04-05-norefit/data001/data001'
datarun2 = load_data(tmp);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2,'load_sta',[],'keep_java_sta',true);
datarun2 = set_polarities(datarun2);
datarun2 = load_ei(datarun2, 'all');




% UDCR
tmp='/Volumes/Analysis/2012-09-24-5/d01-04-05-norefit/data004/data004'
datarun = load_data(tmp);
datarun = load_params(datarun,'verbose',1);
datarun = set_polarities(datarun);
datarun = load_ei(datarun, 'all');
datarun = load_neurons(datarun);

stimulus=read_stim_lisp_output_ath('2011-12-13-2','s07');
parsed=parse_stim_rgbs_ath(stimulus);
map=load(parsed.mappath);
figure
imagesc(map)
% find stimuli patterns
myStim=zeros(4,length(stimulus.rgbs));
for i=1:length(stimulus.rgbs) 
    myStim(1:4,i)=stimulus.rgbs{i}(:,1);
end

[B, ic, ipat] = unique(myStim', 'rows'); % 73 patterns, 40 repetitions


cellTypes=zeros(1,length(datarun2.cell_ids));
for j=1:length(datarun2.cell_types)
    if ~isempty(datarun2.cell_types{j}.cell_ids)
        [~,ia, ib]=intersect(datarun2.cell_ids,datarun2.cell_types{j}.cell_ids);
        cellTypes(ia)=j;
    end
end

for myCell=1:length(datarun.cell_ids)
    % calculate response
    spikes=datarun.spikes{myCell}*1000;
    tr=datarun.triggers(1:2:end)*1000;
            
    myCellType=datarun2.cell_types{cellTypes(myCell)}.name;

    tit=[int2str(myCell),', ID ',int2str(datarun.cell_ids(myCell)), ', ',myCellType];
    
    mySpikes=[];
    for i=1:length(ic)
        patts=find(ipat==i);
        cnt=1;
        for j=1:length(patts)
            mySpikes{cnt,i}=spikes(spikes>tr(patts(j))&spikes<tr(patts(j))+750)-tr(patts(j));
            cnt=cnt+1;
        end
    end
    
    figure
    set(gcf,'position',[1           1        1920        1105])
    for i=1:4
        for j=1:4
            if i==j
                patN=find(B(:,i)==0.48&sum(B,2)==0.48);
                patN2=find(B(:,i)==0.24&sum(B,2)==0.24);
                patN3=find(B(:,i)==-0.48&sum(B,2)==-0.48);
            else
                patN=find(B(:,i)==0.48&B(:,j)==-0.48);
                patN2=find(B(:,i)==0.24&B(:,j)==-0.48);
                patN3=find(B(:,i)==0.12&B(:,j)==-0.48);
            end
            
            
            t=[];
            cnt=0;
            fr=0;
            fr2=0;fr3=0;
            for k=1:40
                t=[t mySpikes{k,patN}'+800*cnt];
                
                tmp=convolved(mySpikes{k,patN},10,800);
                fr=fr+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                tmp=convolved(mySpikes{k,patN2},10,800);
                fr2=fr2+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                tmp=convolved(mySpikes{k,patN3},10,800);
                fr3=fr3+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2));
                
                cnt=cnt+1;
            end
            sbind=sub2ind([4,4],j,i);
            h=subplot(4,4,sbind);
            hold on
            plot(fr/20,'r','linewidth',2)
            plot(fr2/20,'b','linewidth',2)
            plot(fr3/20,'g','linewidth',2)
            if sbind==1
                hl=legend('1','0.5','-1');
                set(hl, 'Box', 'off')
                set(hl, 'Color', 'none')
                set(hl,'location','best')
                title(tit)
            elseif sbind==2
                hl=legend('1,-1','0.5,-1','0.25,-1');
                set(hl, 'Box', 'off')
                set(hl, 'Color', 'none')
                set(hl,'location','best')
            end
            rasterplot(t,cnt,800,h)
            
        end
    end
    
    saveas(gcf,['/Users/alexth/Desktop/ONmidgets_facilitation/',date,'/cell_', ...
        int2str(myCell),'_ID_',int2str(datarun.cell_ids(myCell)),'.jpg'])
    close all
    
end