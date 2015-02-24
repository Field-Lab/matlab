date='2011-12-13-2';

% STA
run='data004-0';
path2load = fullfile(server_path(), [date, '/streamed/',run,'/',run]);
datarunsta = load_data(path2load);
datarunsta = load_params(datarunsta,'verbose',1);
datarunsta = load_sta(datarunsta,'load_sta',[],'keep_java_sta',true);
datarunsta = set_polarities(datarunsta);
datarunsta = load_ei(datarunsta, 'all');


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

a = map_ei(datarun, datarunsta);



cellID=6601;
myCellSTA=find(datarunsta.cell_ids==cellID);
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
        myCellType=['STA ID ',int2str(tmp),', ',datarunsta.cell_types{cellTypes(2,tmp)}.name];
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

cellTypes=zeros(2,length(datarunsta.cell_ids));
cellTypes(1,:)=datarunsta.cell_ids;
for j=1:length(datarunsta.cell_types)
    if ~isempty(datarunsta.cell_types{j}.cell_ids)
        [~,ia, ib]=intersect(cellTypes(1,:),datarunsta.cell_types{j}.cell_ids);
        cellTypes(2,ia)=j;
    end
end

cellMap=zeros(size(datarun.cell_ids));cnt=1;
for myCell=1:length(datarun.cell_ids)
    if isempty(a{myCell})
        cellMap(cnt)=0;
    else
        cellMap(cnt)=a{myCell};
        find(datarunsta.cell_ids==cellMap(cnt));
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
datarunsta = load_data(tmp);
datarunsta = load_params(datarunsta,'verbose',1);
datarunsta = load_sta(datarunsta,'load_sta',[],'keep_java_sta',true);
datarunsta = set_polarities(datarunsta);
datarunsta = load_ei(datarunsta, 'all');


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


cellTypes=zeros(1,length(datarunsta.cell_ids));
for j=1:length(datarunsta.cell_types)
    if ~isempty(datarunsta.cell_types{j}.cell_ids)
        [~,ia, ib]=intersect(datarunsta.cell_ids,datarunsta.cell_types{j}.cell_ids);
        cellTypes(ia)=j;
    end
end

for myCell=1:length(datarun.cell_ids)
    % calculate response
    spikes=datarun.spikes{myCell}*1000;
    tr=datarun.triggers(1:2:end)*1000;
            
    myCellType=datarunsta.cell_types{cellTypes(myCell)}.name;

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
datarunsta = load_data(tmp);
datarunsta = load_params(datarunsta,'verbose',1);
datarunsta = load_sta(datarunsta,'load_sta',[],'keep_java_sta',true);
datarunsta = set_polarities(datarunsta);
datarunsta = load_ei(datarunsta, 'all');


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


cellTypes=zeros(1,length(datarunsta.cell_ids));
for j=1:length(datarunsta.cell_types)
    if ~isempty(datarunsta.cell_types{j}.cell_ids)
        [~,ia, ib]=intersect(datarunsta.cell_ids,datarunsta.cell_types{j}.cell_ids);
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
            
    myCellType=datarunsta.cell_types{cellTypes(myCell)}.name;

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
date='2012-09-13-2';


% STA
run='data001-from-data001_data004_data005';
tmp='/Volumes/Analysis/2012-09-13-2/d01-04-05-norefit/data001-from-data001_data004_data005/data001-from-data001_data004_data005';
datarunsta = load_data(tmp);
datarunsta = load_params(datarunsta,'verbose',1);
datarunsta = load_sta(datarunsta,'load_sta','all','keep_java_sta',true);
datarunsta = set_polarities(datarunsta);
datarunsta = load_neurons(datarunsta);


% UDCR
run='data004-from-data001_data004_data005';
tmp='/Volumes/Analysis/2012-09-13-2/d01-04-05-norefit/data004-from-data001_data004_data005/data004-from-data001_data004_data005';
datarun = load_data(tmp);
datarun = load_neurons(datarun);


% STA 2
run='data005-from-data001_data004_data005';
tmp='/Volumes/Analysis/2012-09-13-2/d01-04-05-norefit/data005-from-data001_data004_data005/data005-from-data001_data004_data005';
datarunsta2 = load_data(tmp);
datarunsta2 = load_params(datarunsta2,'verbose',1);
datarunsta2 = load_sta(datarunsta2,'load_sta','all','keep_java_sta',true);
datarunsta2 = set_polarities(datarunsta2);
datarunsta2 = load_neurons(datarunsta2);




stimulus=read_stim_lisp_output_ath('2012-09-13-2','s04');
parsed=parse_stim_rgbs_ath(stimulus);
map=load('/Volumes/Analysis/2012-09-13-2/stimuli/1234d01/map-0000.txt');
figure
imagesc(map)

myboundaries={};
for i=1:4
    BW=map;
    BW(BW~=i)=0;
    boundaries = bwboundaries(BW);
    myboundaries=[myboundaries boundaries];
end

figure
imagesc(map)
col='rgby'
hold on
for i=1:4
    for k=1:size(myboundaries,1)
        b = myboundaries{k,i};
        plot(b(:,2),b(:,1),col(i));
    end
end

% find stimuli patterns
myStim=zeros(4,length(stimulus.rgbs));
for i=1:length(stimulus.rgbs) 
    myStim(1:4,i)=stimulus.rgbs{i}(:,1);
end

[B, ic, ipat] = unique(myStim', 'rows'); % 73 patterns, 40 repetitions


cellTypes=zeros(1,length(datarunsta.cell_ids));
for j=1:length(datarunsta.cell_types)
    if ~isempty(datarunsta.cell_types{j}.cell_ids)
        [~,ia, ib]=intersect(datarunsta.cell_ids,datarunsta.cell_types{j}.cell_ids);
        cellTypes(ia)=j;
    end
end

    
coord_tform = coordinate_transform(datarunsta,'sta');
ctr=nan(length(datarun.cell_ids),2);
rad=ctr;
fit_angle=nan(length(datarun.cell_ids),1);
for cell_index = 1:length(datarunsta.cell_ids)
    % get the fit
    the_fit = datarunsta.stas.fits{cell_index};
    % skip if doesn't exist
    if isempty(the_fit);continue;end
    % get center
    ctr(cell_index,:) = the_fit.mean;
    % get center radius
    rad(cell_index,:) = the_fit.sd;
    fit_angle(cell_index)=the_fit.angle;
end

sbc=[];
for i=datarunsta.cell_types{8}.cell_ids
    sbc=[sbc find(datarunsta.cell_ids==i)];
end

file_path=['/Users/alexth/Desktop/ONmidgets_facilitation/',date,'/crap/'];
for myCell=sbc
    % calculate response
    spikes=datarun.spikes{myCell}*1000;
    tr=datarun.triggers(1:2:end)*1000;
            
    myCellType=datarunsta.cell_types{cellTypes(myCell)}.name;

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
    
    
    fig=figure('Visible','off','PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[1           1        1920        1105]);
    
    cnt=1;
    for j=0.75:-0.24:0.02
        for i=0.05:0.22:0.7
            h(cnt)=subplot('position',[i j 0.2 0.2]);
            set(h(cnt),'xtick',0,'ytick',25);
            cnt=cnt+1;
        end
    end
    h(13)=subplot('position',[0.75 0.6 0.2 0.3]);
    set(h(13),'xtick',0,'ytick',0);
    h(14)=subplot('position',[0.75 0.2 0.2 0.3]);
    set(h(14),'xtick',0,'ytick',0);

    sbind=1;maxx=50;
    for i=1:4
        for j=setdiff(1:4,i)
            
            pat1pos=find(B(:,i)==0.48&sum(B,2)==0.48); % single cone pos stim
            pat1neg=find(B(:,i)==-0.48&sum(B,2)==-0.48); % single cone neg stim
            pat2=find(B(:,i)==0.48&B(:,j)==-0.48); % pair opposite signs
            
        
            t=[];
            cnt=0;
            fr1pos=0;
            fr1neg=0;fr2=0;
            for k=1:40
                t=[t mySpikes{k,pat2}'+800*cnt];
                
                tmp=convolved(mySpikes{k,pat1pos},10,800);
                fr1pos=fr1pos+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2))/40;
                tmp=convolved(mySpikes{k,pat1neg},10,800);
                fr1neg=fr1neg+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2))/40;
                tmp=convolved(mySpikes{k,pat2},10,800);
                fr2=fr2+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2))/40;
                
                cnt=cnt+1;
            end
            
            maxx=max([fr1pos fr1neg fr2 maxx]);
            
            subplot(h(sbind))            
            hold on
            
            plot(fr1pos,'b','linewidth',2)
            plot(fr1neg,'r','linewidth',2)
            plot(fr2,'g','linewidth',2)

                leg1pos=['+/0'];
                leg1neg=['-/0'];
                leg2=['+/-'];
                hl=legend(leg1pos,leg1neg, leg2);
                set(hl, 'Box', 'off')
                set(hl, 'Color', 'none')
                set(hl,'location','best')
                subtit=['cone ', int2str(i), ' vs cone ', int2str(j)];
                title(subtit)
            
%             title(tit)
%             rasterplot(t,cnt,800,h(sbind))
                        
            sbind=sbind+1;
        end
    end
     

    for sbind=1:12
        subplot(h(sbind))
        axis([0 750 0 maxx])
    end
    
    tmp=double(squeeze(datarunsta.stas.stas{myCell}));
    tmp=imresize(tmp(:,:,5),2, 'method','nearest');
    tmp=tmp/(max(abs(tmp(:)))*2)+0.5;
 
    subplot(h(13))
    colormap gray
    imagesc(tmp);
    hold on

    for m=1:4
        for k=1:size(myboundaries,1)
            b = myboundaries{k,m};
            plot(b(:,2),b(:,1),col(m));
        end
    end

 
    [X, Y] = drawEllipse([ctr(myCell,:)*2 rad(myCell,:)*4 fit_angle(myCell)]);
    [X, Y] = tformfwd(coord_tform, X, Y);
%     [y,x]=find(map>0);
%     
%     IN = inpolygon(x,y,X,Y);
%     figure
%     plot(X,Y,'y')
%     hold on
%     imagesc(map)
%     find(IN)
%     x(find(IN))
%     y(find(IN))
%     figure
%     imagesc(map())

    minx=max([min(X) 0]);
    maxx=min([max(X) 600]);
    miny=max([min(Y) 0]);
    maxy=min([max(Y) 600]);
    if maxx<0; maxx=100; end
           
    axis([minx maxx miny maxy])

      
    tmp=double(squeeze(datarunsta2.stas.stas{myCell}));
    tmp=imresize(tmp(:,:,5),2, 'method','nearest');
    tmp=tmp/(max(abs(tmp(:)))*2)+0.5;
 
    subplot(h(14))
    colormap gray
    imagesc(tmp);
    hold on
    for m=1:4
        for k=1:size(myboundaries,1)
            b = myboundaries{k,m};
            plot(b(:,2),b(:,1),col(m));
        end
    end
    axis([minx maxx miny maxy])
    
    
    tit=['cell_',int2str(myCell),'_ID_',int2str(datarun.cell_ids(myCell))];
    
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',file_path,tit));
    
    close(fig)
    
end


%% 2012-09-24-5 d01-04-05
date='2012-09-24-5';

% STA
tmp='/Volumes/Analysis/2012-09-24-5/d01-04-05-norefit/data001/data001'
datarunsta = load_data(tmp);
datarunsta = load_params(datarunsta,'verbose',1);
datarunsta = load_sta(datarunsta,'load_sta','all','keep_java_sta',true);
datarunsta = set_polarities(datarunsta);
datarunsta = load_ei(datarunsta, 'all');


% UDCR
tmp='/Volumes/Analysis/2012-09-24-5/d01-04-05-norefit/data004/data004';
datarun = load_data(tmp);
datarun = load_neurons(datarun);


% STA 2
tmp='/Volumes/Analysis/2012-09-24-5/d01-04-05-norefit/data005/data005';
datarunsta2 = load_data(tmp);
datarunsta2 = load_params(datarunsta2,'verbose',1);
datarunsta2 = load_sta(datarunsta2,'load_sta','all','keep_java_sta',true);
datarunsta2 = set_polarities(datarunsta2);
datarunsta2 = load_neurons(datarunsta2);




stimulus=read_stim_lisp_output_ath('2012-09-24-5','s04');
parsed=parse_stim_rgbs_ath(stimulus);
map=load('/Volumes/Analysis/2012-09-24-5/stimuli/1234d01/map-0000.txt');
figure
imagesc(map)

myboundaries={};
for i=1:4
    BW=map;
    BW(BW~=i)=0;
    boundaries = bwboundaries(BW);
    myboundaries=[myboundaries boundaries];
end

figure
imagesc(map)
col='rgby'
hold on
for i=1:4
    for k=1:size(myboundaries,1)
        b = myboundaries{k,i};
        plot(b(:,2),b(:,1),col(i));
    end
end

% find stimuli patterns
myStim=zeros(4,length(stimulus.rgbs));
for i=1:length(stimulus.rgbs) 
    myStim(1:4,i)=stimulus.rgbs{i}(:,1);
end

[B, ic, ipat] = unique(myStim', 'rows'); % 73 patterns, 40 repetitions


cellTypes=zeros(1,length(datarunsta.cell_ids));
for j=1:length(datarunsta.cell_types)
    if ~isempty(datarunsta.cell_types{j}.cell_ids)
        [~,ia, ib]=intersect(datarunsta.cell_ids,datarunsta.cell_types{j}.cell_ids);
        cellTypes(ia)=j;
    end
end

    
coord_tform = coordinate_transform(datarunsta,'sta');
ctr=nan(length(datarun.cell_ids),2);
rad=ctr;
fit_angle=nan(length(datarun.cell_ids),1);
for cell_index = 1:length(datarunsta.cell_ids)
    % get the fit
    the_fit = datarunsta.stas.fits{cell_index};
    % skip if doesn't exist
    if isempty(the_fit);continue;end
    % get center
    ctr(cell_index,:) = the_fit.mean;
    % get center radius
    rad(cell_index,:) = the_fit.sd;
    fit_angle(cell_index)=the_fit.angle;
end

sbc=[];
for i=datarunsta.cell_types{3}.cell_ids
    sbc=[sbc find(datarunsta.cell_ids==i)];
end

file_path=['/Users/alexth/Desktop/ONmidgets_facilitation/',date,'/on_midget/'];
for myCell=sbc
    % calculate response
    spikes=datarun.spikes{myCell}*1000;
    tr=datarun.triggers(1:2:end)*1000;
            
    myCellType=datarunsta.cell_types{cellTypes(myCell)}.name;

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
    
    
    fig=figure('Visible','off','PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[1           1        1920        1105]);
    
    cnt=1;
    for j=0.75:-0.24:0.02
        for i=0.05:0.22:0.7
            h(cnt)=subplot('position',[i j 0.2 0.2]);
            set(h(cnt),'xtick',0,'ytick',25);
            cnt=cnt+1;
        end
    end
    h(13)=subplot('position',[0.75 0.6 0.2 0.3]);
    set(h(13),'xtick',0,'ytick',0);
    h(14)=subplot('position',[0.75 0.2 0.2 0.3]);
    set(h(14),'xtick',0,'ytick',0);

    sbind=1;maxx=50;
    for i=1:4
        for j=setdiff(1:4,i)
            
            pat1pos=find(B(:,i)==0.48&sum(B,2)==0.48); % single cone pos stim
            pat1neg=find(B(:,i)==-0.48&sum(B,2)==-0.48); % single cone neg stim
            pat2=find(B(:,i)==0.48&B(:,j)==-0.48); % pair opposite signs
            
        
            t=[];
            cnt=0;
            fr1pos=0;
            fr1neg=0;fr2=0;
            for k=1:40
                t=[t mySpikes{k,pat2}'+800*cnt];
                
                tmp=convolved(mySpikes{k,pat1pos},10,800);
                fr1pos=fr1pos+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2))/40;
                tmp=convolved(mySpikes{k,pat1neg},10,800);
                fr1neg=fr1neg+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2))/40;
                tmp=convolved(mySpikes{k,pat2},10,800);
                fr2=fr2+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2))/40;
                
                cnt=cnt+1;
            end
            
            maxx=max([fr1pos fr1neg fr2 maxx]);
            
            subplot(h(sbind))            
            hold on
            
            plot(fr1pos,'b','linewidth',2)
            plot(fr1neg,'r','linewidth',2)
            plot(fr2,'g','linewidth',2)

                leg1pos=['+/0'];
                leg1neg=['-/0'];
                leg2=['+/-'];
                hl=legend(leg1pos,leg1neg, leg2);
                set(hl, 'Box', 'off')
                set(hl, 'Color', 'none')
                set(hl,'location','best')
                subtit=['cone ', int2str(i), ' vs cone ', int2str(j)];
                title(subtit)
            
%             title(tit)
%             rasterplot(t,cnt,800,h(sbind))
                        
            sbind=sbind+1;
        end
    end
     
    
    for sbind=1:12
        subplot(h(sbind))
        axis([0 750 0 maxx])
    end
    
    tmp=double(squeeze(datarunsta.stas.stas{myCell}));
    tmp=imresize(tmp(:,:,5),2, 'method','nearest');
    tmp=tmp/(max(abs(tmp(:)))*2)+0.5;
 
    subplot(h(13))
    colormap gray
    imagesc(tmp);
    hold on

    for m=1:4
        for k=1:size(myboundaries,1)
            b = myboundaries{k,m};
            plot(b(:,2),b(:,1),col(m));
        end
    end

 
    [X, Y] = drawEllipse([ctr(myCell,:)*2 rad(myCell,:)*4 fit_angle(myCell)]);
    [X, Y] = tformfwd(coord_tform, X, Y);
%     [y,x]=find(map>0);
%     
%     IN = inpolygon(x,y,X,Y);
%     figure
%     plot(X,Y,'y')
%     hold on
%     imagesc(map)
%     find(IN)
%     x(find(IN))
%     y(find(IN))
%     figure
%     imagesc(map())

    minx=max([min(X) 0]);
    maxx=min([max(X) 600]);
    miny=max([min(Y) 0]);
    maxy=min([max(Y) 600]);
    if maxx<0; maxx=100; end
           
    axis([minx maxx miny maxy])

      
    tmp=double(squeeze(datarunsta2.stas.stas{myCell}));
    tmp=imresize(tmp(:,:,5),2, 'method','nearest');
    tmp=tmp/(max(abs(tmp(:)))*2)+0.5;
 
    subplot(h(14))
    colormap gray
    imagesc(tmp);
    hold on
    for m=1:4
        for k=1:size(myboundaries,1)
            b = myboundaries{k,m};
            plot(b(:,2),b(:,1),col(m));
        end
    end
    axis([minx maxx miny maxy])
    
    
    tit=['cell_',int2str(myCell),'_ID_',int2str(datarun.cell_ids(myCell))];
    
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',file_path,tit));
    
    close(fig)
    
end




%% 2012-09-24-5 d03-06-07
date='2012-09-24-5';
path2data='/Volumes/Analysis/2012-09-24-5/d03-06-07-norefit/'
starun='data003';
sta2run='data007';
udrun='data006';
stim='s06';
stimmap='/Volumes/Analysis/2012-09-24-5/stimuli/1234d03/map-0000.txt';


date='2012-09-13-2';
path2data='/Volumes/Analysis/2012-09-13-2/d01-04-05-norefit/';
starun='data001-from-data001_data004_data005';
sta2run='data005-from-data001_data004_data005';
udrun='data004-from-data001_data004_data005';
stim='s04';
stimmap='/Volumes/Analysis/2012-09-13-2/stimuli/1234d01/map-0000.txt';



date='2011-12-13-2';
path2data='/Volumes/Analysis/2011-12-13-2/data004_data007_data008-norefit/';
starun='data004-from-data004_data007_data008';
sta2run='data008-from-data004_data007_data008';
udrun='data007-from-data004_data007_data008';
stim='s07';
stimmap='/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f04_1234/map-0000.txt';



% STA
datarunsta = load_data(fullfile(path2data,starun,starun));
datarunsta = load_params(datarunsta,'verbose',1);
datarunsta = load_sta(datarunsta,'load_sta','all','keep_java_sta',true);
datarunsta = set_polarities(datarunsta);

% STA 2
datarunsta2 = load_data(fullfile(path2data,sta2run,sta2run));
datarunsta2 = load_params(datarunsta2,'verbose',1);
datarunsta2 = load_sta(datarunsta2,'load_sta','all','keep_java_sta',true);
datarunsta2 = set_polarities(datarunsta2);

% UDCR
datarun = load_data(fullfile(path2data,udrun,udrun));
datarun = load_neurons(datarun);


% read stimulus info
stimulus=read_stim_lisp_output_ath(date,stim);
parsed=parse_stim_rgbs_ath(stimulus);
map=load(stimmap);

myboundaries=cell(14,4);
for i=1:4
    BW=map;
    BW(BW~=i)=0;
    boundaries = bwboundaries(BW);
    myboundaries(1:size(boundaries,1),i)=boundaries;
end
figure
imagesc(map)
col='rgby';
hold on
for i=1:4
    for k=1:size(myboundaries,1)
        b = myboundaries{k,i};
        if ~isempty(b)
            plot(b(:,2),b(:,1),col(i));
        end
    end
end

% find stimuli patterns
myStim=zeros(4,length(stimulus.rgbs));
for i=1:length(stimulus.rgbs) 
    myStim(1:4,i)=stimulus.rgbs{i}(:,1);
end
[B, ic, ipat] = unique(myStim', 'rows'); % 73 patterns, 40 repetitions

% find cell types
cellTypes=zeros(1,length(datarunsta.cell_ids));
for j=1:length(datarunsta.cell_types)
    if ~isempty(datarunsta.cell_types{j}.cell_ids)
        [~,ia, ib]=intersect(datarunsta.cell_ids,datarunsta.cell_types{j}.cell_ids);
        cellTypes(ia)=j;
    end
end
cellTypeList=unique(cellTypes);

% 
coord_tform = coordinate_transform(datarunsta,'sta');
ctr=nan(length(datarun.cell_ids),2);
rad=ctr;
fit_angle=nan(length(datarun.cell_ids),1);
for cell_index = 1:length(datarunsta.cell_ids)
    % get the fit
    the_fit = datarunsta.stas.fits{cell_index};
    % skip if doesn't exist
    if isempty(the_fit);continue;end
    % get center
    ctr(cell_index,:) = the_fit.mean;
    % get center radius
    rad(cell_index,:) = the_fit.sd;
    fit_angle(cell_index)=the_fit.angle;
end


curType='rest';

switch curType
    case 'sbc'
        allCells=find(cellTypes==5);
        pol=-1;
    case 'on_midget'
        allCells=find(cellTypes==3);
        pol=1;
    case 'off_midget'
        allCells=find(cellTypes==4);
        pol=-1;
    case 'on_parasol'
        allCells=find(cellTypes==1);
        pol=1;
    case 'off_parasol'
        allCells=find(cellTypes==2);
        pol=-1;
    case 'rest'
        allCells=find(cellTypes>5);
        pol=1;
end

file_path=['/Users/alexth/Desktop/ONmidgets_facilitation/',date,'/',curType,'/'];
if ~exist(file_path,'dir')
    mkdir(file_path);
end

leg1pos='+/0';
leg1neg='-/0';
fullSTA=1;
for myCell=allCells
    % calculate response
    spikes=datarun.spikes{myCell}*1000;
    tr=datarun.triggers(1:2:end)*1000;
            
    myCellType=datarunsta.cell_types{cellTypes(myCell)}.name;

    tit=[int2str(myCell),', ID ',int2str(datarun.cell_ids(myCell)), ', ',myCellType];
    
    mySpikes=[];
    for i=1:length(ic)
        patts=find(ipat==i);
        cnt=1;
        for j=1:length(patts)
            if length(tr)>= patts(j)
                mySpikes{cnt,i}=spikes(spikes>tr(patts(j))&spikes<tr(patts(j))+750)-tr(patts(j));                
            else
                mySpikes{cnt,i}=[];
            end
            cnt=cnt+1;
        end
    end
    
    
    fig=figure('visible','off','PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[1           1        1920        1105]);
    
    cnt=1;
    for j=0.75:-0.24:0.02
        for i=0.05:0.22:0.7
            h(cnt)=subplot('position',[i j 0.2 0.2]);
            set(h(cnt),'xtick',0,'ytick',25);
            cnt=cnt+1;
        end
    end
    h(13)=subplot('position',[0.75 0.6 0.2 0.3]);
    set(h(13),'xtick',0,'ytick',0);
    h(14)=subplot('position',[0.75 0.2 0.2 0.3]);
    set(h(14),'xtick',0,'ytick',0);

    sbind=1;maxx=50;
    for i=1:4
        for j=setdiff(1:4,i)
            
            pat1pos=find(B(:,i)==0.48&sum(B,2)==0.48); % single cone pos stim
            pat1neg=find(B(:,i)==-0.48&sum(B,2)==-0.48); % single cone neg stim
            if pol==-1
                pat2=find(B(:,i)==-0.48&B(:,j)==0.48); % pair opposite signs
                leg2='-/+';
            else
                pat2=find(B(:,i)==0.48&B(:,j)==-0.48); % pair opposite signs
                leg2='+/-';
            end
            
        
            t=[];
            cnt=0;
            fr1pos=0;
            fr1neg=0;fr2=0;
            for k=1:40
                t=[t mySpikes{k,pat2}'+800*cnt];
                
                tmp=convolved(mySpikes{k,pat1pos},10,800);
                fr1pos=fr1pos+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2))/40;
                tmp=convolved(mySpikes{k,pat1neg},10,800);
                fr1neg=fr1neg+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2))/40;
                tmp=convolved(mySpikes{k,pat2},10,800);
                fr2=fr2+tmp(((size(tmp,2)-800)/2+1):end-((size(tmp,2)-800)/2))/40;
                
                cnt=cnt+1;
            end
            
            maxx=max([fr1pos fr1neg fr2 maxx]);
            
            subplot(h(sbind))            
            hold on
            
            plot(fr1pos,'b','linewidth',2)
            plot(fr1neg,'r','linewidth',2)
            plot(fr2,'g','linewidth',2)


                hl=legend(leg1pos,leg1neg, leg2);
                set(hl, 'Box', 'off')
                set(hl, 'Color', 'none')
                set(hl,'location','best')
                subtit=['cone ', int2str(i), ' vs cone ', int2str(j)];
                title(subtit)
            
%             title(tit)
%             rasterplot(t,cnt,800,h(sbind))
                        
            sbind=sbind+1;
        end
    end
     
    
    for sbind=1:12
        subplot(h(sbind))
        axis([0 750 0 maxx])
    end
    
    tmp=double(squeeze(datarunsta.stas.stas{myCell}));
    tmp=imresize(tmp(:,:,5),2, 'method','nearest');
    tmp=tmp/(max(abs(tmp(:)))*2)+0.5;
 
    subplot(h(13))
    colormap gray
    imagesc(tmp);
    hold on

    for m=1:4
        for k=1:size(myboundaries,1)
            b = myboundaries{k,m};
            if ~isempty(b)
                plot(b(:,2),b(:,1),col(m));
            end
        end
    end

 
    [X, Y] = drawEllipse([ctr(myCell,:)*2 rad(myCell,:)*4 fit_angle(myCell)]);
    [X, Y] = tformfwd(coord_tform, X, Y);
%     [y,x]=find(map>0);
%     
%     IN = inpolygon(x,y,X,Y);
%     figure
%     plot(X,Y,'y')
%     hold on
%     imagesc(map)
%     find(IN)
%     x(find(IN))
%     y(find(IN))
%     figure
%     imagesc(map())

    minx=max([min(X) 0]);
    maxx=min([max(X) 600]);
    miny=max([min(Y) 0]);
    maxy=min([max(Y) 600]);
    if maxx<0; maxx=100; end
    if maxy<0; maxy=100; end
           
    if fullSTA
        axis([0 600 0 600])
    else
        axis([minx maxx miny maxy])
    end

      
    tmp=double(squeeze(datarunsta2.stas.stas{myCell}));
    tmp=imresize(tmp(:,:,5),2, 'method','nearest');
    tmp=tmp/(max(abs(tmp(:)))*2)+0.5;
 
    subplot(h(14))
    colormap gray
    imagesc(tmp);
    hold on
    for m=1:4
        for k=1:size(myboundaries,1)
            b = myboundaries{k,m};
            if ~isempty(b)
                plot(b(:,2),b(:,1),col(m));
            end
        end
    end
    if fullSTA
        axis([0 600 0 600])
    else
        axis([minx maxx miny maxy])
    end
    
    
    tit=['cell_',int2str(myCell),'_ID_',int2str(datarun.cell_ids(myCell))];
    
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',file_path,tit));
    
    close(fig)
    
end



