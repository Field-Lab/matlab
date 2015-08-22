%% 2011-12-13-2
date='2011-12-13-2';


% STA
run='data000-from-d00_04';
tmp='/Volumes/Analysis/2011-12-13-2/d00-04-norefit/data000-from-d00_04/data000-from-d00_04';
datarunsta = load_data(tmp);
datarunsta = load_params(datarunsta,'verbose',1);
datarunsta = load_sta(datarunsta,'load_sta','all','keep_java_sta',true);
datarunsta = set_polarities(datarunsta);
datarunsta = load_neurons(datarunsta);


% UDCR
run='data003-from-d00_04';
tmp='/Volumes/Analysis/2011-12-13-2/d00-04-norefit/data003-from-d00_04/data003-from-d00_04';
datarun = load_data(tmp);
datarun = load_neurons(datarun);


% STA 2
run='data004-from-d00_04';
tmp='/Volumes/Analysis/2011-12-13-2/d00-04-norefit/data004-from-d00_04/data004-from-d00_04';
datarunsta2 = load_data(tmp);
datarunsta2 = load_params(datarunsta2,'verbose',1);
datarunsta2 = load_sta(datarunsta2,'load_sta','all','keep_java_sta',true);
datarunsta2 = set_polarities(datarunsta2);
datarunsta2 = load_neurons(datarunsta2);




stimulus=read_stim_lisp_output_ath('2011-12-13-2','s03');
parsed=parse_stim_rgbs_ath(stimulus);
map=load('/Volumes/Data/2011-12-13-2/Visual/2011-12-13-2_f00_1234/map-0000.txt');
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

offm=[];
for i=datarunsta.cell_types{4}.cell_ids
    offm=[offm find(datarunsta.cell_ids==i)];
end

file_path=['/Users/alexth/Desktop/ONmidgets_facilitation/',date,'_data003/off_midget/'];
if ~isdir(file_path)
    mkdir(file_path)
end

for myCell=offm
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
