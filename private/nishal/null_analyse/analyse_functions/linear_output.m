function gen=linear_output(WN_datafile,InterestingCell_vis_id,ref_cell_number,condMovies,sta_to_use,nConditions,cond_str)


%% Load STA
opt=struct('load_all',true);
datarun=load_data(WN_datafile,opt)
datarun=load_sta(datarun)
datarun=load_params(datarun)
datarun=load_ei(datarun,'all','array_type',519);
%get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
%STAs!

vision_id=InterestingCell_vis_id(ref_cell_number);
idx=[1:length(datarun.cell_ids)];
matlab_id=idx(datarun.cell_ids==vision_id);
cell_ei=datarun.ei.eis{matlab_id};
ssta=datarun.stas.stas{matlab_id};

%% Fit STA

if(strcmp(sta_to_use,'fit'))
    
addpath(genpath('~/GITs/matlab/code')); 
addpath(genpath('~/GITs/matlab/private/nishal/fwdfittingfunctions'));
fit_info = fit_sta(ssta)
%datarun = compute_sta_fits(datarun, 'all');

%params=datarun.matlab.sta_fits{matlab_id}.initial_params;
params=fit_info.initial_params;
full_fit = sta_fit_function(params);

figure;
for itime=30:-1:1
    itime
    subplot(2,1,1);
    imagesc(squeeze(sum(ssta(:,:,:,itime),3)));
    colormap gray
    axis image
    colorbar
    caxis([min(ssta(:)), max(ssta(:))]);

    
    subplot(2,1,2);
    imagesc(squeeze(sum(full_fit(:,:,:,itime),3)));
    colormap gray
    axis image
    colorbar 
    caxis([min(ssta(:)), max(ssta(:))]);
 
     pause(1/120);
    
end
end


%%


if(strcmp(sta_to_use,'fit'))
s_use=full_fit(:,:,:,end:-1:1);
s_use(:,:,:,15:end)=0; % Clip fitted STA also!
end

if(strcmp(sta_to_use,'orig'))
s_use=ssta(:,:,:,end:-1:1);
end

if(strcmp(sta_to_use,'clip'))
    stas{1}=ssta(:,:,:,end:-1:1);
    cell_params.STAlen=14;
    stas_use = clipSTAs(stas,cell_params);
    s_use=stas_use{1};
end

%s_use(:,:,:,15:end)=0; % DOUBT .. TODO ? Ask EJ
%%


gen=cell(nConditions,1);
for icond=1:nConditions
% Format movie into right format
mov_toformat=condMovies{icond};
mov=zeros(size(mov_toformat,2),size(mov_toformat,3),3,size(mov_toformat,1));

for itime=1:size(mov,4)
mov(:,:,1,itime)=squeeze(mov_toformat(itime,:,:));
mov(:,:,2,itime)=squeeze(mov_toformat(itime,:,:));
mov(:,:,3,itime)=squeeze(mov_toformat(itime,:,:));
end
color_l='rkbmcg';
% Generate response

movie_time=size(mov,4);

gen_signals=zeros(movie_time,3);
for col=1:3;
    col
    st_temp=zeros(size(s_use,2),size(s_use,1),1,size(s_use,4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:30
        st_temp(:,:,:,itime)=s_use(:,:,col,itime)'; % DOUBT .. Could be a reason for things to fail!!!!!
    end
    s_use_new=st_temp;

Filtlen = size(s_use_new,4);
Filtdim1=size(s_use_new,1);
Filtdim2=size(s_use_new,2);

movie_new_len=movie_time;
mov2=zeros(Filtdim1 ,Filtdim2,movie_new_len+Filtlen-1);
mov2(:,:,Filtlen:movie_new_len+Filtlen-1)=squeeze(mov(:,:,col,:)); % Append zeros before the movie
sz=max(size(mov2,3)-size(s_use_new,4) + 1, 0);

gen_signals(:,col) = reshape(convn(mov2,squeeze(s_use_new(end:-1:1,end:-1:1,1,:)),'valid'),[sz,1]);
end
    %
gen{icond}=sum(gen_signals,2); 
end


figure;
for icond=1:nConditions
subplot(nConditions,1,icond);
plot([1:length(gen{icond})]/120,gen{icond},color_l(icond));
ylim([min(gen{1}) , max(gen{1})]);
title(sprintf(' Vision ID: %d Condition : %s STA type %s',InterestingCell_vis_id(ref_cell_number),cond_str{icond},sta_to_use));

end
end