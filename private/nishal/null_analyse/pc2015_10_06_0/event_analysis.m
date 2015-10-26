%% Sub-unit estimation 

addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%%
% Condition strings
nConditions=6;
condDuration=1272/6;
cond_str=cell(6,1);
cond_str{1}='Original';
cond_str{2}='On selected ';
cond_str{3}='All ON';
cond_str{4}='Off selected ';
cond_str{5}='All OFF';
cond_str{6}='SBC';

interestingConditions=[1,2,3,4,5,6];
%% Load Movies


% make movies
interval=2;
condMov=cell(nConditions,1);
rawMovFrames=1272/(2);
icnt=0;
% make pixel histogram
for imov=[1,2,4,6,8,10]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-2/Visual/pc2015_03_09_2_data038/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    qq=permute(movie,[2,3,1]);
    ifcnt = 0;
    condMov{icnt}=zeros(size(qq,1),size(qq,2),size(qq,3)*interval);
    for iframe=1:size(qq,3)
        for irepeat=1:interval
            ifcnt=ifcnt+1;
            condMov{icnt}(:,:,ifcnt)=qq(:,:,iframe); % cond mov is between 0 and 1 now!
        end
        
    end
    
end

%%
WN_datafile = '2015-03-09-2/streamed/data038/data038';
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-2/data038/data038';
neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];


datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[2]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

% InterestingCell_vis_id = [4548]%[662,6826,4548] %NullCells1;
condDuration=10.6;
nConditions=6;
GLM_fit_link= '/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_ON parasol/';
ref_cell_number=22;
    cellID=InterestingCell_vis_id(ref_cell_number);
    
    % make directory
       if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number))))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data042/CellType_%s/CellID_%d/',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)));
       end
       
    % Plot recorded raster
     [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);
       
   event_tol=0.03;
   thr=0.4;
    [spkCondColl,eventOverlap,h_event] = event_count(spkCondColl,condDuration,thr,event_tol); 
    InterestingCell_vis_id(ref_cell_number)
    %  pause

%% have spkCondColl, spike rate during spiking event of originial conditions

ref_cond=1;
null_cond=5;
null_event_rate=[];
ref_event_rate=spkCondColl(ref_cond).events_spks;
bin_select=zeros(size(spkCondColl(ref_cond).rec_rast,2),1);
null_event_rate_bins = bin_select;

for ievent=1:length(spkCondColl(ref_cond).events_list_start)
bin_start = spkCondColl(ref_cond).events_bin_start(ievent);
bin_end = spkCondColl(ref_cond).events_bin_end(ievent);
bin_select(bin_start:bin_end)=1;
null_event_rate(ievent)=sum(sum(spkCondColl(null_cond).rec_rast(:,bin_start:bin_end)))/(size(spkCondColl(null_cond).rec_rast,1)*spkCondColl(ref_cond).events_list_length(ievent));
null_event_rate_bins(bin_start:bin_end) = null_event_rate(ievent);
end
figure;
hist(null_event_rate,100);
xlim([0,80])

%% 
idx=1:length(datarun.cell_ids);
matlab_id=idx(datarun.cell_ids==cellID);
time_c = datarun.vision.timecourses(matlab_id).r +  datarun.vision.timecourses(matlab_id).g +  datarun.vision.timecourses(matlab_id).b;
time_c=(time_c(end:-1:1));
filtMov=cell(nConditions,1);
for icond=1:nConditions
filtMov{icond}=zeros(size(condMov{icond},1),size(condMov{icond},2),size(condMov{icond},3));
    for x=1:size(condMov{icond},1)
        for y=1:size(condMov{icond},2)
        filtMov{icond}(x,y,:)=conv(squeeze(condMov{icond}(x,y,:)),time_c,'same');
        end
    end
end


%%
datarun=load_sta(datarun);
stas{1}=double(datarun.stas.stas{matlab_id});


stas_new=cell(length(stas),1);
for icell=1:length(stas)
    st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:size(stas{1},4)
        st_temp(:,:,:,itime)=mean(stas{icell}(:,:,:,end-itime+1),3)'; % DOUBT .. Could be a reason for things to fail!!!!!
    end
    stas_new{icell}=st_temp;
end
stas=stas_new;

% Used in movie post process
cell_params.STAlen=14;
[stas_clipped,totalMaskAccept,CellMasks]= clipSTAs(stas,cell_params);

%% Filter movie using mask
maskedMov=cell(nConditions,1);

for icond=1:nConditions
maskedMov{icond}=zeros(sum(CellMasks{1}(:)),size(filtMov{icond},3));

    for itime=1:size(filtMov{icond},3)
     xx=filtMov{icond}(:,:,itime);
     maskedMov{icond}(:,itime) = xx(logical(CellMasks{1}));
    end
end

%% weights
icond=5;
rec_rast=spkCondColl(icond).rec_rast;
binSz=10;
w=zeros(size(rec_rast,2)/binSz,1);
bin_select_coarse=zeros(size(rec_rast,2)/binSz,1);
event_rate_coarse=zeros(size(rec_rast,2)/binSz,1);
for itime=1:size(rec_rast,2)
iframe = floor((itime-1)/binSz)+1;
w(iframe) = w(iframe)+sum(rec_rast(:,itime));
bin_select_coarse(iframe) = bin_select_coarse(iframe) + bin_select(itime);
event_rate_coarse(iframe) = null_event_rate_bins(itime);
end

bin_select_coarse=bin_select_coarse>0;
%% Optimization problem
k=3;
icond=1;
stim_dim=size(maskedMov{icond},1);
X = maskedMov{icond};

X=X(:,bin_select_coarse);
w=event_rate_coarse(bin_select_coarse);

th=1;
lam=0;
%coeffs= 2*(double(w-th > 0)-0.5); 
coeffs = w-th;

%%
cvx_begin
variable a(stim_dim,k)
minimize (coeffs'*X'*a(:,1) + coeffs'*X'*a(:,2) +coeffs'*X'*a(:,3)+lam*norm(a(:),1)) 
sum(a,2)<=ones(stim_dim,1)
cvx_end

%% Spectral Clustering on spike rate! 

npts=size(X,2);
A=size(npts,npts);
sigma=0.35;
%sigma=5;
for i=1:npts
    for j=1:npts
  %  A(i,j) = exp(-(w(i)-w(j))^2/(2*sigma^2));
   A(i,j) = exp(-norm(X(:,i)-X(:,j))^2/(2*sigma^2));
    end
end

%
d=sum(A,2);
D=diag(d);

L=eye(npts) - (D^(-0.5))*A*(D^(-0.5));

[E,DE]=eig(L);
lam=diag(DE);
[lam_sort,I]=sort(lam,'ascend');

E_sort=E(:,I);

figure;
plot(abs(lam_sort),'*');
title('Eigen value spectrum');

k=3;
Y = E_sort(:,1:k);
U=0*Y;
for i=1:npts
U(i,:) = Y(i,:)/norm(Y(i,:));
end

[idx,C]=kmeans(U,k);

figure;
scatter(U(:,1),U(:,2),[],idx+3,'filled');
title('K-means clustering');

% Cluster centers in original space.

Centers=zeros(size(X,1),k);
for ik=1:k
Centers(:,ik)=mean(X(:,idx==ik),2);
end

figure;
plot(Centers)