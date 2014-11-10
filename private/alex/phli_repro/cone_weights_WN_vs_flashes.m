%% load datarun with white noise
piece = '2012-09-13-2';
run = 'data009';

% define data path
datarunA = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunA=load_data(datarunA,opt);
datarunA = load_params(datarunA,struct('verbose',1));  
datarunA = load_sta(datarunA,'load_sta',[]);
datarunA = set_polarities(datarunA);
datarunA = load_cones(datarunA,'bayes');
datarunA = make_mosaic_struct(datarunA);
datarunA = get_sta_fits_from_vision(datarunA);  
datarunA = make_voronoi_masks(datarunA);


%% load datarun with single flashes

piece = '2012-09-13-2';
run = 'data013';

datarun = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt = struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarun = load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  

offP_d09 = [376 1471 4636 4771 5551];
offM_d09 = [151 1396 1741 1816 1863 1966 2026 2116 6796 6991];

%% map cells in dataruns (single flashes as master)

[cell_list, failed_cells]=map_ei(datarun,datarunA);

my_cells=zeros(length(cell_list),2);
for i=1:length(cell_list)
    my_cells(i,1)=datarun.cell_ids(i);
    if isempty(cell_list{i})
        my_cells(i,2)=0;
    else
        my_cells(i,2)=cell_list{i};
    end
end



[mapped_cells_offP_d09, islave, imaster]=intersect(offP_d09,my_cells(:,2));

[mapped_cells_offM_d09, islave, imaster]=intersect(offM_d09,my_cells(:,2));

%% load stimulus and maps

stimulus=read_stim_lisp_output_ath('2012-09-13-2','s13');

parsed=parse_stim_rgbs_ath(stimulus);

load(['/Volumes/Analysis/2012-09-13-2/stimuli/allconesd09/map-0000.txt'])
figure
imagesc(map_0000)
length(unique(map_0000))
a=parsed.rgbs;


k=cell(1,stimulus.numcones); % stimuli list for each cone
m=zeros(size(stimulus.pulses,2)/stimulus.numcones,stimulus.numcones); % when stimuli were applied
for i=1:size(stimulus.pulses,2)
    p=a{i};
    t=find(p(:,1));
    k{t}=[k{t} p(t,1)];
    p=find(m(:,t)==0,1);
    m(p,t)=i;
end

% m - n_of_indices x repetition_number

%% plot single flashes responses

tr=datarun.triggers(1:2:end);
mmin=0;mmax=0.75;

myMaps=zeros(600,600,length(imaster));
cnt=1;
for i=imaster'
    i
    kmaps=map_0000;
    spikes=datarun.spikes{i};
    conv=zeros(29,750);
    maxFiringRate=zeros(29,2);
    maxFiringRateLP=zeros(29,2);
    figure
    for j=1:29
        trig_local=m(k{j}==-0.48,j);
        for f=trig_local'
            h=spikes-tr(f);
            hh= h>=mmin & h<mmax;
            hh=round(h(hh)*1000);
            tmp=convolved(hh,40,750);
            tmp=tmp(121:end-120);
            conv(j,:)=conv(j,:)+tmp/50;
        end
        subplot(5,6,j)
        plot(conv(j,:))
        conv(j,:)=conv(j,:)-mean(conv(j,500:600));
        [maxFiringRate(j,1) maxFiringRateLP(j,1)]=max(conv(j,:));
        hold on
        
        conv(j,:)=conv(j,:)-conv(j,:);
        trig_local=m(k{j}==-0.288,j);
        for f=trig_local'
            h=spikes-tr(f);
            hh= h>=mmin & h<mmax;
            hh=round(h(hh)*1000);
            tmp=convolved(hh,40,750);
            tmp=tmp(121:end-120);
            conv(j,:)=conv(j,:)+tmp/50;
        end
        plot(conv(j,:),'r')
        conv(j,:)=conv(j,:)-mean(conv(j,500:600));
        [maxFiringRate(j,2) maxFiringRateLP(j,2)]=max(conv(j,:));
        
        kmaps(map_0000==j)=sum(maxFiringRate(j,:));
        axis([0 750 0 30])
    end
    subplot(5,6,30)    
    plot(maxFiringRate)
    
    kmaps=kmaps/max(kmaps(:));
    myMaps(:,:,cnt)=kmaps;
    figure
    imagesc(kmaps)
    hold on
    plot_rf_summaries(datarunA, mapped_cells_offM_d09(cnt), 'scale',2,'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')
    
    cnt=cnt+1;

  
end


%% plot WN responses and cone maps

kmaps=map_0000;
cnt=1;
for i=imaster'
    i
    myCell=mapped_cells_offM_d09(cnt);
    
    datarunA = load_sta(datarunA,'load_sta',myCell);
    
    sta=-double(squeeze(datarunA.stas.stas{datarunA.cell_ids==myCell}));
    sta=sta(:,:,5);
    sta=imresize(sta,2);
    
    sta=sta/max(abs(sta(:)))/2;
    sta=sta+0.5;
    
    kmaps(kmaps>0)=1;
     
    figure
    subplot(1,2,1)
    colormap gray
    imagesc(sta) 
    
    merged=zeros(600,600,3);
    merged(:,:,2)=sta;
    merged(:,:,1)=kmaps/2;
    
    subplot(1,2,2)
    imagesc(merged)
 
    cnt=cnt+1;
  
end

%% get cone map for stimulation

masks = datarunA.cones.mosaic.voronoi_masks;
excludes=[];
indexes = 1:length(masks);
indexes = setdiff(indexes, excludes);

min_neighbor_dist = 2;
max_self_dist     = 3.5;
spaced = space_voronoi_masks(datarunA, min_neighbor_dist, max_self_dist);
cone_map = index_masks(spaced, num2cell(indexes));

imagesc(cone_map)


%% plot cone weights
cnt=1;
for i=imaster'
    i
    myCell=mapped_cells_offM_d09(cnt);
    
    datarunA = load_sta(datarunA,'load_sta',myCell);
    
    coneWeights=datarunA.cones.weights(:,datarunA.cell_ids==myCell);
    coneWeights=coneWeights/max(coneWeights);
    [sortedWeights,coneID]=sort(coneWeights,'descend');
    coneID=coneID(1:10);
    
    tmp=zeros(600,600);
    for j=1:10
        tmp(cone_map==coneID(j))=sortedWeights(j);
    end
       
    tmp=-tmp/2;
    
    myConeResponse=myMaps(:,:,cnt)/2;
    
    
    to_plot=tmp+myConeResponse+0.5;
    to_plot(1,1)=0;
    to_plot(1,2)=1;
    
    figure
    colormap gray
    imagesc(to_plot)
    hold on
    plot_rf_summaries(datarunA, mapped_cells_offM_d09(cnt), 'scale',2,'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')

    
%     merged=zeros(600,600,3);
%     merged(:,:,2)=tmp;
%     merged(:,:,1)=myConeResponse;
%  
%     figure
%     imagesc(merged) 
    cnt=cnt+1;
  
end

