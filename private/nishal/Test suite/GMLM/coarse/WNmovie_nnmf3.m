%% GMLM with MEL and EM for normal cells at coarser resolution .. 

location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library

addpath(('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act2'));
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/GLM/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/GLM2/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/code'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act_2/'));
%% Dataset details
% WN_datafile = '2015-03-09-2/data031/data031';
% WN_datafile_short='2015-03-09-2/data031/data031';
% movie_xml = 'BW-8-6-0.48-11111-40x40';
% stim_length=1800;% in seconds

WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
stim_length=1800;% 
%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);
%%
for cellID = [datarun.cell_types{12}.cell_ids,datarun.cell_types{13}.cell_ids,datarun.cell_types{14}.cell_ids,datarun.cell_types{5}.cell_ids];
  
    cellID


%cell_glm_fit = sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_OFF parasol/CellID_%d.mat',cellID);
%load(cell_glm_fit);
%% Set-up parameters
extract_movie_response2;

 for nSU=[5,3,2,4,7,6,8,9]
     
%load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/CellID_%d/Cell%d_full_su4_5.mat',cellID,cellID));
 %% nnmf
 mov_use=maskedMovdd;
 binnedResponsesbigd = spksGen;

 STX = mov_use(:,binnedResponsesbigd~=0);
 
 STX = STX-min(STX(:));
 STX = STX/max(STX(:));
 
     su_log=zeros(size(STX,1),size(STX,1));
     
     W_log=cell(50,1);
     H_log=cell(50,1);
 for ifit=1:50
     close all
     ifit
     [W,H] = nnmf(STX,nSU);
 
 W_log{ifit}=W;
 H_log{ifit}=H;
 %% plot
h= figure('Color','w');
   mask = totalMaskAccept;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

u_spatial_log = zeros(40,40,nSU);

for ifilt=1:nSU

u_spatial = reshape_vector(W(:,ifilt)/sum(W(:,ifilt)),masked_frame,indexedframe);
%u_spatial = reshape_vector(WN_uSq{ifilt},masked_frame,indexedframe);

subplot(ceil((nSU+1)/2),ceil((nSU+1)/2),ifilt)
imagesc(u_spatial(x_coord,y_coord));
colormap gray
colorbar
title(sprintf('nnmf Filter: %d',ifilt));
axis square
u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;
end

% [v,I] = max(u_spatial_log,[],3);
% 
% xx=I.*(v>0.2);
% xx=xx(x_coord,y_coord);
% 
% subplot(2,2,nSU);
% imagesc(repelem(xx,10,10));
% title(sprintf('Num SU : %d',nSU));
% hold on;

if ~exist(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Large Cells/CellID_%d/nnmf/',cellID),'dir');
    mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Large Cells/CellID_%d/nnmf/',cellID));
end

if ~exist(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Large Cells/CellID_%d/nnmf/SU_%d',cellID,nSU),'dir');
    mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Large Cells/CellID_%d/nnmf/SU_%d',cellID,nSU));
end


hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Large Cells/CellID_%d/nnmf/SU_%d/Cell%d_nnmf_su_%d_fit_%d.eps',cellID,nSU,cellID,nSU,ifit));

%% make occurance histogram
 k_est=W';
 su=double((k_est)==repmat(max((k_est),[],1),[nSU,1]));
 su_log = su_log + su'*su;
 end
 %% make near-neighbour co-occurance map
 h=figure;
 ilist= 1:length(masked_frame);
 for ix=x_coord
    for iy=y_coord
    i1 =ilist(indexedframe(ix,iy)==masked_frame);
    if(~isempty(i1))
   
    i2 = ilist(indexedframe(ix,iy+1)==masked_frame);
    
    if(~isempty(i2))
        if(su_log(i1,i2)>0)
    plot([ix,ix],[iy,iy+1],'LineWidth',su_log(i1,i2));
    hold on;
     text(ix,iy+0.5,sprintf('%d',su_log(i1,i2)));
    hold on;
        end
    end
    
    i2 = ilist(indexedframe(ix,iy-1)==masked_frame);
    if(~isempty(i2))
          if(su_log(i1,i2)>0)
    plot([ix,ix],[iy,iy-1],'LineWidth',su_log(i1,i2));
    hold on;
       text(ix,iy-0.5,sprintf('%d',su_log(i1,i2)));
    hold on;
          end
    end
    
    i2 = ilist(indexedframe(ix-1,iy)==masked_frame);
    if(~isempty(i2))
           if(su_log(i1,i2)>0)
    plot([ix-1,ix],[iy,iy],'LineWidth',su_log(i1,i2));
    hold on;
       text(ix-0.5,iy,sprintf('%d',su_log(i1,i2)));
    hold on;
           end
    end
    
    
    i2 = ilist(indexedframe(ix+1,iy)==masked_frame);
    if(~isempty(i2))
           if(su_log(i1,i2)>0)
    plot([ix+1,ix],[iy,iy],'LineWidth',su_log(i1,i2));
    hold on;
       text(ix+0.5,iy,sprintf('%d',su_log(i1,i2)));
    hold on;
           end
    end
    
    end
    
    end
 end
 
 hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Large Cells/CellID_%d/nnmf/Cell%d_nnmf_su_%d.eps',cellID,cellID,nSU));
save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Large Cells/CellID_%d/nnmf/Cell%d_nnmf_su_%d.mat',cellID,cellID,nSU),'su_log','W_log','H_log','STX');

 end
 end
 