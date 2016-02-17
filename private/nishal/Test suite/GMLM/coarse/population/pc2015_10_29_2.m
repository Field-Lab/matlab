
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library

addpath(('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act2'));
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('/Volumes/cd Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/GLM/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/GLM2/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/code'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act_2/'));
%% Dataset details

WN_datafile = '2015-10-29-2/d00_36-norefit/data002/data002'; % data002 RGB-8-2-0.24-11111, 30 min, all streamed , low contrast (0.24)
WN_datafile_short=WN_datafile;
movie_xml = 'RGB-8-2-0.48-11111';
stim_length=1800;% 
%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

icell_l=0;
cells = [datarun.cell_types{1}.cell_ids];
ttf_log = zeros(30,length(cells));
totalMaskAccept_log = zeros(80*40,length(cells));
Y = zeros(length(cells),216000);

for cellID = cells
    icell_l=icell_l+1
extract_movie_response4;
ttf_log(:,icell_l) =ttf;
totalMaskAccept_log(:,icell_l) = totalMaskAccept(:);
Y(icell_l,:) = spksGen;
end

 ttf_avg = mean(ttf_log,2);

 mov=squeeze(mean(mov,3));
 maskedMovdd= filterMov(mov,ones(size(mov,1),size(mov,2)),squeeze(ttf_avg));

  save('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_analysis/On_parasol.mat','maskedMovdd','Y','ttf_log','ttf_avg','totalMaskAccept_log','cells','-v7.3')
  
  
  %% fit ASM
  
path = '/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_analysis'
load([path,'/Off_parasol.mat']);
binnedSpikeResponses_coll = Y;
total_mask_log = totalMaskAccept_log;

cells = double(cells);
cellsChoose = (cells ==3287) | (cells ==3318 ) | (cells ==3155) | (cells ==3066);
cellsChoose=logical(cellsChoose);
mask = sum(total_mask_log(:,cellsChoose),2)~=0;
ttf = ttf_avg;


ifit=0;
for fitNnum=1:1
for Ns=[4,8,6,10,12,14,16];
%fitASM_pop = fitASM_EM_Population(1000*maskedMovdd,binnedSpikeResponses_coll(cellsChoose,:),Ns,mask);
ifit=ifit+1
lam_start = 0.1;
gamma1=0.0000;
gamma2 = 0;
initVal=[];
%initVal.K = maskedMovdd(mask,:)*binnedSpikeResponses_coll(cellsChoose,:)';
%initVal.B = diag(ones(Ns,1));
[fitASM_pop,fval] = fitASM_EM_Population_sparse_split_admm(maskedMovdd,binnedSpikeResponses_coll(cellsChoose,:),Ns,mask,gamma1,gamma2,lam_start,initVal);
B_use= plotSU_withcells(fitASM_pop.K,mask,total_mask_log(:,cellsChoose),exp(fitASM_pop.B));
save(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_%d.mat',ifit),'Ns','fitASM_pop','fval','mask','gamma1','gamma2','lam_start','initVal','cellsChoose');
pause(0.5);
end
end 
  
%% Do decoding

%% Load generic data

path = '/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_analysis'
load([path,'/Off_parasol.mat']);
binnedSpikeResponses_coll = Y;
total_mask_log = totalMaskAccept_log;


cells = double(cells);
cellsChoose = (cells ==3287) | (cells ==3318 ) | (cells ==3155) | (cells ==3066);
cellsChoose=logical(cellsChoose);
mask = sum(total_mask_log(:,cellsChoose),2)~=0;
ttf = ttf_avg;

% calculate STAs
Y = binnedSpikeResponses_coll;
ncell_type1 = size(Y,1);

stas_celltype1 = cell(ncell_type1,1);

for icell=1:ncell_type1
stas_celltype1{icell} = maskedMovdd*Y(icell,:)'/sum(Y(icell,:));
end



% Load fits
for ifit=2%1:7
    ifit
load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_%d.mat',ifit));
B_use= plotSU_withcells(fitASM_pop.K,reshape(mask,[80,40]),total_mask_log(:,cellsChoose),exp(fitASM_pop.B));

% Do decoding 
T=4000;
rho=1000; % for rho=1, it worked ok ? 
lambda=0.1;
Xref = maskedMovdd(mask,1:T);

%[Xdecode,errMap] = decode_ASM_population3D(fitASM_pop,binnedSpikeResponses_coll(cellsChoose,1:100),mask,ttf,rho,lambda,1000*maskedMovdd(mask,1:100),B_use)
spks = binnedSpikeResponses_coll(cellsChoose,1:T);
[Xdecode_ifit,errMap] = decode_ASM_population(fitASM_pop,spks,mask,ttf,rho,lambda,Xref,B_use)
save(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_%d_decode_experimental.mat',ifit),'spks','ifit','Xref','Xdecode_ifit','cellsChoose','T','ttf','B_use','mask','rho','lambda','fitASM_pop')
end

%% Analyze decoding

%% load decodings
dec1 = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_1_decode.mat');
dec2 = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_2_decode.mat');
dec3 = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_3_decode.mat');
celltype_sign = -1; % -1 for off and +1 for on
load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_analysis/Off_stas.mat')

% 
 path = '/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_analysis'
 load([path,'/Off_parasol.mat'],'totalMaskAccept_log');

 ifit=3
 load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_%d.mat',ifit));
B_use= plotSU_withcells(fitASM_pop.K,reshape(mask,[80,40]),totalMaskAccept_log(:,cellsChoose),exp(fitASM_pop.B));

data_frndly = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_1_2_3_decode_nonGPU.mat');
B_use= plotSU_withcells(data_frndly.K3,reshape(mask,[80,40]),totalMaskAccept_log(:,cellsChoose),exp(data_frndly.B3));

[B_use,h,B_use_lr]= plotSU_withcells(data_frndly.K1,reshape(mask,[80,40]),totalMaskAccept_log(:,cellsChoose),exp(data_frndly.B1));

%% T kernel
T = 4000;
ttf = dec1.ttf;
Tker =sparse(T,T);
for itime=1:T
for jtime = max(itime-29,1):itime
    Tker(itime,jtime) = ttf((itime-jtime)+1)/norm(ttf);
end

% make it circulant matrix
if(itime<30)
Tker(itime,end-30+itime:end) = ttf(end:-1:end-30+itime)/norm(ttf);
end
end
%Tkerinv = (pinv(0.01*eye(T,T)+ Tker));

%% 
X_coll{1}  = dec1.Xref;
X_coll{2} = dec1.Xdecode_ifit;
X_coll{3} = dec2.Xdecode_ifit;
X_coll{4} = dec3.Xdecode_ifit;
nMov = length(X_coll);
mask = dec1.mask;
cell_sig=[1,1,1,1];

mov_dec=cell(nMov);
s1=80;s2=40;
for imov=1:length(X_coll)
mov_show = gather((Tker\((cell_sig(imov)*X_coll{imov})'))'); 
mov_show = temporal_filter(mov_show,ttf,0.7);
[mov_show,xlim1,ylim1]= movie_reshape(mov_show,s1,s2,mask);
% mov_show1 = mov_show1(xlim1(1):xlim1(2),ylim1(1):ylim1(2),:);
mov_dec{imov} = spatial_filter(mov_show,stas_celltype1,0.25,1,s1,s2);
end

xx= zeros(s1,s2);xx(xlim1(1):xlim1(2),ylim1(1):ylim1(2))=1;xx=logical(xx);
mask_mov = repelem(xx,1,1,T);

figure;
for itime=2200:2210
   for imov=1:nMov
   subplot(1,nMov,imov);
   imagesc(mov_dec{imov}(:,:,itime)');axis image;colormap gray;colorbar;%caxis([min(mov_dec_reshape(:)), 0.5*max(mov_dec_reshape(:))]);
   ylim(xlim1);xlim(ylim1);
   end
   pause
   %pause(1/120);
end

tidx=2400:2800;
figure;jpix = 12; ipix=59;
plot(squeeze(mov_dec{1}(ipix,jpix,tidx)));
hold on;
for imov=2:nMov
scale=norm(squeeze(mov_dec{1}(ipix,jpix,tidx)))/norm(squeeze(mov_dec{imov}(ipix,jpix,tidx)));
plot(squeeze(scale*mov_dec{imov}(ipix,jpix,tidx)));
end

% find correlations
%% Find correlation between pixel time courses
tidx = 1000:3000;
mov2 = mov_dec{1};
corr_log=cell(nMov,1);
for imov=1:nMov
mov1 = mov_dec{imov};

corr_pix=zeros(size(mov1,1),size(mov1,2));

for idim1=1:size(mov1,1)
    for idim2 = 1:size(mov1,2)
        seq1 = squeeze(mov1(idim1,idim2,tidx));
        
        seq2 = squeeze(mov2(idim1,idim2,tidx));
        corr_pix(idim1,idim2) = corr(seq1,seq2);
    end
end
corr_log{imov}=(corr_pix');
end

figure;
icnt=0;
for imov=1:nMov
    for jmov=1:nMov
    icnt=icnt+1;
        subplot(nMov,nMov,icnt);
        if(imov==jmov)
        imagesc((corr_log{imov}));axis image;ylim(xlim1);xlim(ylim1);colormap gray
        else
        imagesc((corr_log{imov}-corr_log{jmov}));axis image;ylim(xlim1);xlim(ylim1);colormap gray 
        end
        title(sprintf('%d-%d',imov,jmov));
    end
end

col='rgbm';
figure;
icnt=0;
for imov=1:nMov
    for jmov=1:nMov
    icnt=icnt+1;
        subplot(nMov,nMov,icnt);
        if(imov==jmov)
        imagesc((corr_log{imov}));axis image;ylim(xlim1);xlim(ylim1);colormap gray
        else
        imagesc(repelem(sqrt(sum((mov_dec{imov}(:,:,tidx)-mov_dec{jmov}(:,:,tidx)).^2,3))',20,20));axis image;ylim(20*xlim1);xlim(20*ylim1);colormap gray 
        
        for i=1:4
        hold on;
        plot(B_use{i}(:,1),B_use{i}(:,2),col(i));
        end
       
        
        end
        title(sprintf('%d-%d',imov,jmov));
    end
end

% maximally differentiated
col='rgbm';
figure;
icnt=0;
iidx=[1:4000]';
for imov=2:nMov
    imov
    for jmov=2:nMov
    icnt=icnt+1;
        subplot(nMov-1,nMov-1,icnt);
        if(imov==jmov)
       
        else
            % find correlations
            
            mov2 = mov_dec{1};
            %for imov
            mov1 = mov_dec{imov};
            corr_pix_imov=zeros(size(mov1,1),size(mov1,2));corr_pix_jmov=corr_pix_imov;
            for idim1=1:size(mov1,1)
                for idim2 = 1:size(mov1,2)
                    
                    % find times
                    diff=squeeze(sum(sum((mov_dec{imov}(idim1,idim2,:)-mov_dec{jmov}(idim1,idim2,:)).^2,1),2));
                    thr = prctile((diff),70);
                    tdiff=(diff>thr)&(iidx>1000&iidx<3000);
                    
                    seqimov = squeeze(mov_dec{imov}(idim1,idim2,tdiff));
                    seq2 = squeeze(mov2(idim1,idim2,tdiff));
                    corr_pix_imov(idim1,idim2) = corr(seqimov,seq2);
                    
                    
                    seqjmov = squeeze(mov_dec{jmov}(idim1,idim2,tdiff));
                    seq2 = squeeze(mov2(idim1,idim2,tdiff));
                    corr_pix_jmov(idim1,idim2) = corr(seqjmov,seq2);
                end
            end
            corr_log1=(corr_pix_imov');
            corr_log2=(corr_pix_jmov');
          
            imagesc(repelem(corr_log1-corr_log2,20,20));axis image;ylim(20*xlim1);xlim(20*ylim1);colormap gray 
        
        
        for i=1:4
        hold on;
        plot(B_use{i}(:,1),B_use{i}(:,2),col(i));
        end
       
        
        end
        title(sprintf('%d - %d',imov,jmov));
    end
end


%% Decode null stimulus

path = '/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_analysis'
load([path,'/Off_parasol.mat']);
binnedSpikeResponses_coll = Y;
total_mask_log = totalMaskAccept_log;


cells = double(cells);
cellsChoose = (cells ==3287) | (cells ==3318 ) | (cells ==3155) | (cells ==3066);
cellsChoose=logical(cellsChoose);
mask = sum(total_mask_log(:,cellsChoose),2)~=0;
ttf = ttf_avg;

% calculate STAs
Y = binnedSpikeResponses_coll;
ncell_type1 = size(Y,1);

stas_celltype1 = cell(ncell_type1,1);

for icell=1:ncell_type1
stas_celltype1{icell} = maskedMovdd*Y(icell,:)'/sum(Y(icell,:));
end

% load null stimulus response and stimulus
% load basic data
nConditions=8;
condDuration=1200/4;
cond_str=cell(nConditions,1);
interestingConditions=[1,2,3,4,5,6,7];

dataRuns_OFF_additivity = [3,4,6,7,9,11,13];
dataRuns_ON_additivity = [3,5,6,8,10,12,13];
movies_OFF_addivitiy =[1,2,5,6,10,14,13];
movies_ON_additivity = [1,4,5,8,12,16,13];
mov_id=6;idata=7;
location = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit';

% make movies
interval=4;
condMov=cell(nConditions,1);
rawMovFrames=1200/(4);
icnt=0;

% make pixel histogram
imov=mov_id
    %[stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-10-29-2/Visual/%d.rawMovie',imov),rawMovFrames,1);
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_from_01/%d.rawMovie',imov),rawMovFrames,1);
    
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
            condMov{icnt}(:,:,ifcnt)=qq(:,:,iframe)+0.5; % cond mov is between 0 and 1 now!
        end
        
    end
    
mov=condMov{icnt};
movd = (mov-0.5);  % *2 if it is a low contrast movie!
maskedMov= filterMov(movd,mask,ttf);
Xref = maskedMov;

% load spikes
condDuration=10;
nConditions=1;
cond_str={};
for itrial =1:30
spks=zeros(4,1200);icell=0;
for cellID = [3287,3318 ,3155,3066]
    icell=icell+1;
Null_datafile = sprintf('%s/data0%02d',location,idata);
neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',idata)];
[spkColl,spkCondColl,h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
spks_cell=makeSpikeMat(spkCondColl.spksColl,1/120,1200);
spks(icell,:)=spks_cell(itrial,:);
end

% Load fits
for ifit=1:7%1:7
    ifit
load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_%d.mat',ifit));
B_use= plotSU_withcells(fitASM_pop.K,reshape(mask,[80,40]),total_mask_log(:,cellsChoose),exp(fitASM_pop.B));

% Do decoding 
T=4000;
rho=1000; % for rho=1, it worked ok ? 
lambda=0.1;
%Xref = maskedMovdd(mask,1:T);
%spks = binnedSpikeResponses_coll(cellsChoose,1:T);

[Xdecode_ifit,errMap] = decode_ASM_population(fitASM_pop,spks,mask,ttf,rho,lambda,Xref,B_use)
save(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_%d_decode_trial_%d_null.mat',ifit,itrial),'spks','ifit','Xref','Xdecode_ifit','cellsChoose','T','ttf','B_use','mask','rho','lambda','fitASM_pop')
end
end

%% Analyze multiple decodings!
%% load decodings
for ifit=1:7
    for itrial=1:30
dec{ifit,itrial} = load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_%d_decode_trial_%d_null.mat',ifit,itrial));
    end
end

celltype_sign = -1; % -1 for off and +1 for on
load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_analysis/Off_stas.mat')

% 
path = '/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_analysis'
load([path,'/Off_parasol.mat'],'totalMaskAccept_log');

data_frndly = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_sharedSU/Off_parasol_1_2_3_decode_nonGPU.mat');
[B_use,h,B_use_lr]= plotSU_withcells(data_frndly.K1,reshape(mask,[80,40]),totalMaskAccept_log(:,cellsChoose),exp(data_frndly.B1));

%% T kernel
T = 1200;
ttf = dec{1,1}.ttf;
Tker =sparse(T,T);
for itime=1:T
for jtime = max(itime-29,1):itime
    Tker(itime,jtime) = ttf((itime-jtime)+1)/norm(ttf);
end

% make it circulant matrix
if(itime<30)
Tker(itime,end-30+itime:end) = ttf(end:-1:end-30+itime)/norm(ttf);
end
end
%Tkerinv = (pinv(0.01*eye(T,T)+ Tker));

%% see all 30 movies and calculate average ?
mov_dec=cell(7,30);s1=80;s2=40;
mean_mov=cell(7,1);
for ifit=1:7
    ifit
    mean_mov{ifit}=zeros(s1,s2,1200);
for itrial=1:30
    X_dec = dec{ifit,itrial}.Xdecode_ifit;
    
mov_show = gather((Tker\((X_dec)'))'); 
mov_show = temporal_filter(mov_show,ttf,0.7);
[mov_show,xlim1,ylim1]= movie_reshape(mov_show,s1,s2,mask);
% mov_show1 = mov_show1(xlim1(1):xlim1(2),ylim1(1):ylim1(2),:);
mov_dec{ifit,itrial} = spatial_filter(mov_show,stas_celltype1,0.25,1,s1,s2);
mean_mov{ifit} = mean_mov{ifit} + mov_dec{ifit,itrial};
close all
end
mean_mov{ifit}=mean_mov{ifit}/30;
end


xx= zeros(s1,s2);xx(xlim1(1):xlim1(2),ylim1(1):ylim1(2))=1;xx=logical(xx);
mask_mov = repelem(xx,1,1,T);

figure;
for itime=500:600
   for itrial=1:30
   subplot(6,5,itrial);
   imagesc(mov_dec{itrial}(:,:,itime)');axis image;colormap gray;colorbar;%caxis([min(mov_dec_reshape(:)), 0.5*max(mov_dec_reshape(:))]);
   ylim(xlim1);xlim(ylim1);
   end
   pause
   %pause(1/120);
end

figure;
for itime=500:600
    for ifit=1:7
        
        subplot(3,3,ifit);
        imagesc(mean_mov{ifit}(:,:,itime)');axis image;colormap gray;colorbar;%caxis([min(mov_dec_reshape(:)), 0.5*max(mov_dec_reshape(:))]);
        ylim(xlim1);xlim(ylim1);
     
    end
    pause;
end

tidx = 200:1000;
icnt=0;
figure;
for ifit=1:7
for jfit =1:7
icnt=icnt+1;
    subplot(7,7,icnt);
    imagesc(sqrt(mean((mean_mov{ifit}(:,:,tidx)-mean_mov{jfit}(:,:,tidx)).^2,3)'));
    axis image;colormap gray;%colorbar;%caxis([min(mov_dec_reshape(:)), 0.5*max(mov_dec_reshape(:))]);
        ylim(xlim1);xlim(ylim1);
        
end
end

% 
x=56;y=17;
figure;
for ifit =[1:7]
plot(tidx,squeeze(mean_mov{ifit}(x,y,tidx)));
hold on;
end


