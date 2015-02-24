date='2014-06-04-7';

% STA
run='data002';
path2load = fullfile(server_path(), [date, '/streamed/',run,'/',run]);
datarun2 = load_data(path2load);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2,'load_sta',[],'keep_java_sta',true);
datarun2 = set_polarities(datarun2);
datarun2 = load_ei(datarun2, 'all');


% UDCR
run='data003';
path2load = fullfile(server_path(), [date, '/',run,'/',run]);
datarun = load_data(path2load);
datarun = load_params(datarun,'verbose',1);
datarun = set_polarities(datarun);
datarun = load_ei(datarun, 'all');
datarun = load_neurons(datarun);

stimulus=read_stim_lisp_output_ath('2014-06-04-7','s03','map-0000');
parsed=parse_stim_rgbs_ath(stimulus);
map=load(parsed.mappath);
figure
imagesc(map)
% find stimuli patterns
myStim=zeros(121,length(stimulus.rgbs));
for i=1:length(stimulus.rgbs) 
    myStim(:,i)=stimulus.rgbs{i}(:,1);
end

[B, ic, ia] = unique(myStim', 'rows'); % 73 patterns, 40 repetitions

a = map_ei(datarun, datarun2);

% find SBC
sbcs2=datarun2.cell_types{1, 5}.cell_ids;
sbcs=[];
for i=1:length(sbcs2)
    for j=1:length(a)
        if a{j}==sbcs2(i)
            sbcs=[sbcs j];
            break;
        elseif j==length(a)
             sbcs=[sbcs 0];             
        end
    end
end

datarun2=load_sta(datarun2,'load_sta',sbcs2);
sta=datarun2.stas.stas{datarun2.cell_ids==sbcs2(1)};
sta=sta(:,:,:,5);
sta=sta/(2*max(sta(:)))+0.5;
tmp=sta(:,:,3);
tmp=tmp-mean(tmp(:));
tmp(tmp<0)=0;
tmp=tmp/max(tmp(:));
mySTA=zeros(320,320,3);
mySTA(:,:,2)=tmp;
tmp=map;
tmp=tmp/max(tmp(:));
% tmp(tmp>0)=1;
mySTA(:,:,1)=tmp;
figure
imagesc(mySTA)




file_path=['/Users/alexth/Desktop/SBC_alex/SBC_single_cone_response/',date,'/'];

for myCell=sbcs
    if myCell>0
        % calculate response
        spikes=datarun.spikes{myCell}*1000;
        tr=datarun.triggers(1:2:end)*1000;
        
        mySpikes=cell(1,size(myStim,2));
        for i=1:size(myStim,2)
            mySpikes{i}=spikes(spikes>tr(i)&spikes<tr(i)+500)-tr(i);
        end
        
        fig=figure('Visible','off','PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
        set(gcf,'color','white','position',[1           1        1920        1105]);

        for i=1:121
            patN=find(myStim(i,:)==0.48);
            
            t=[];
            fr=0;
            for k=1:25
                t=[t mySpikes{patN(k)}'+500*(k-1)];
                
                tmp=convolved(mySpikes{patN(k)},10,500);
                fr=fr+tmp(((size(tmp,2)-500)/2+1):end-((size(tmp,2)-500)/2));
            end
            h=subplot(11,11,i);
            hold on
            plot(fr/20,'r','linewidth',2)
            rasterplot(t,25,500,h)
        end
        
        cnt=1;
        for i=[108 102 95 100 94]
            patN=find(myStim(i,:)==0.48);
            
            t=[];
            fr=0;
            for k=1:25
                t=[t mySpikes{patN(k)}'+500*(k-1)];
                
                tmp=convolved(mySpikes{patN(k)},10,500);
                fr=fr+tmp(((size(tmp,2)-500)/2+1):end-((size(tmp,2)-500)/2));
            end
            h=subplot(2,3,cnt);
            hold on
            plot(fr/20,'r','linewidth',2)
            rasterplot(t,25,500,h)
            cnt=cnt+1;
        end
        
        
        
        tit=['cell_',int2str(datarun.cell_ids(myCell)),'_ID_',int2str(a{myCell})];
        
        print(fig,'-dpdf',sprintf('%s%s%s.pdf',file_path,tit));
        close all
    end
    
end




date='2014-06-04-7';

% STA
run='data002';
path2load = fullfile(server_path(), [date, '/',run,'/',run]);
datarun2 = load_data(path2load);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2,'load_sta',[],'keep_java_sta',true);
datarun2 = set_polarities(datarun2);
datarun2 = load_ei(datarun2, 'all');


run='data004';
path2load = fullfile(server_path(), [date, '/',run,'/',run]);
datarun = load_data(path2load);
datarun = load_params(datarun,'verbose',1);
datarun = set_polarities(datarun);
datarun = load_ei(datarun, 'all');
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);



a = map_ei(datarun, datarun2);


% find SBC
sbcs2=datarun2.cell_types{1, 5}.cell_ids;
sbcs=[];
for i=1:length(sbcs2)
    for j=1:length(a)
        if a{j}==sbcs2(i)
            sbcs=[sbcs j];
            break;
        elseif j==length(a)
             sbcs=[sbcs 0];             
        end
    end
end


datarun2=load_sta(datarun2,'load_sta',sbcs2);

datarun=load_sta(datarun,'load_sta',datarun.cell_ids(sbcs(sbcs>0)));

for mc=1
    
    sta=datarun2.stas.stas{datarun2.cell_ids==sbcs2(mc)};
    sta=sta(:,:,:,5);
    sta=sta/(2*max(sta(:)))+0.5;
    tmp=sta(:,:,3);
    tmp=tmp-mean(tmp(:));
    tmp(tmp<0)=0;
    tmp=tmp/max(tmp(:));
    
    mySTA=zeros(320,320,3);
    mySTA(:,:,1)=tmp;
    
    sta=datarun.stas.stas{sbcs(mc)};
    sta=sta(:,:,:,5);
    sta=sta/(2*max(sta(:)))+0.5;
    tmp=sta(:,:,3);
    tmp=tmp-mean(tmp(:));
    tmp(tmp<0)=0;
    tmp=tmp/max(tmp(:));
    
    mySTA(:,:,3)=tmp;
    
    tmp=map;
    % tmp=tmp/max(tmp(:));
    tmp(tmp>0)=0.3;
    mySTA(:,:,2)=tmp;
    
    figure
    imagesc(mySTA)
end
