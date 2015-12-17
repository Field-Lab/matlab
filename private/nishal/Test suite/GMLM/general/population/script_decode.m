
%% Load generic data
path = '/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/'
load([path,'/Off_type1.mat']);
binnedSpikeResponses_coll = Y;
total_mask_log = totalMaskAccept_log;

cellsChoose = zeros(size(binnedSpikeResponses_coll,1),1);
cell_choose_num = [3,7,8];
cellsChoose(cell_choose_num)=1; % [25,26,31,23] or [25]
cellsChoose=logical(cellsChoose);
mask = sum(total_mask_log(:,cellsChoose),2)~=0;
ttf = ttf_avg;

%% calculate STAs
Y = binnedSpikeResponses_coll;
ncell_type1 = size(Y,1);

stas_celltype1 = cell(ncell_type1,1);

for icell=1:ncell_type1
stas_celltype1{icell} = maskedMovdd*Y(icell,:)'/sum(Y(icell,:));
end



%% Load fits
ifit=17;
load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/Off_type1_fit_%d.mat',ifit));
B_use= plotSU_withcells(fitASM_pop.K,mask,total_mask_log(:,cellsChoose),exp(fitASM_pop.B));

%% Do decoding 
T=4000;
rho=10; % for rho=1, it worked ok ? 
lambda=0.1;
Xref = maskedMovdd(mask,1:T);

%[Xdecode,errMap] = decode_ASM_population3D(fitASM_pop,binnedSpikeResponses_coll(cellsChoose,1:100),mask,ttf,rho,lambda,1000*maskedMovdd(mask,1:100),B_use)
[Xdecode_ifit17,errMap] = decode_ASM_population(fitASM_pop,binnedSpikeResponses_coll(cellsChoose,1:T),mask,ttf,rho,lambda,Xref,B_use)
 

%% Load fits
ifit=1;
load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/Off_type1_fit_%d.mat',ifit));
B_use= plotSU_withcells(fitASM_pop.K,mask,total_mask_log(:,cellsChoose),exp(fitASM_pop.B));

%% Do decoding 

T=4000;
rho=10; % for rho=1, it worked ok ? 
lambda=0.1;
Xref = maskedMovdd(mask,1:T);

%[Xdecode,errMap] = decode_ASM_population3D(fitASM_pop,binnedSpikeResponses_coll(cellsChoose,1:100),mask,ttf,rho,lambda,1000*maskedMovdd(mask,1:100),B_use)
[Xdecode_ifit1,errMap] = decode_ASM_population(fitASM_pop,binnedSpikeResponses_coll(cellsChoose,1:T),mask,ttf,rho,lambda,Xref,B_use)
 
%% T kernel

Tker =gpuArray(zeros(T,T));
for itime=1:T
for jtime = max(itime-29,1):itime
    Tker(itime,jtime) = ttf((itime-jtime)+1)/norm(ttf);
end

% make it circulant matrix
if(itime<30)
Tker(itime,end-30+itime:end) = ttf(end:-1:end-30+itime)/norm(ttf);
end
end
Tkerinv = (pinv(0.01*eye(T,T)+ Tker));
%% Analyze decoding
s1=40;s2=40;

mov_show = gather((Tkerinv*(Xdecode)')'); 
mov_show = temporal_filter(mov_show,ttf,0.7);
[mov_show,xlim1,ylim1]= movie_reshape(mov_show,s1,s2,mask);
% mov_show1 = mov_show1(xlim1(1):xlim1(2),ylim1(1):ylim1(2),:);
mov_dec_reshape = spatial_filter(mov_show,stas_celltype1,0.25,1,s1,s2);

s1=40;s2=40;
mov_show = gather((Tkerinv*Xref')'); 
mov_show = temporal_filter(mov_show,ttf,0.7);
mov_show= movie_reshape(mov_show,s1,s2,mask);
mov_ori_reshape = spatial_filter(mov_show,stas_celltype1,0.25,1,s1,s2);

xx= zeros(s1,s2);xx(xlim1(1):xlim1(2),ylim1(1):ylim1(2))=1;xx=logical(xx);
mask_mov = repelem(xx,1,1,T);

figure;
for itime=200:210%size(Xref,2)
   subplot(1,2,1);
   imagesc(mov_dec_reshape(:,:,itime));axis image;colormap gray;colorbar;%caxis([min(mov_dec_reshape(:)), 0.5*max(mov_dec_reshape(:))]);
   ylim(xlim1);xlim(ylim1);
   
   subplot(1,2,2);
   imagesc(mov_ori_reshape(:,:,itime));axis image;colormap gray;colorbar;%caxis([min(mov_ori_reshape(:)),0.5*max(mov_ori_reshape(:))]);
   ylim(xlim1);xlim(ylim1);
    
   %pause
   pause(1/120);
end

tidx=400:800;
figure;ipix = 12; jpix=20;
plot(squeeze(mov_ori_reshape(ipix,jpix,tidx)));
hold on;
scale=norm(squeeze(mov_ori_reshape(ipix,jpix,tidx)))/norm(squeeze(mov_dec_reshape(ipix,jpix,tidx)));
plot(squeeze(scale*mov_dec_reshape(ipix,jpix,tidx)));


% find correlations
%% Find correlation between pixel time courses
tidx = 300:700;
mov1 = mov_dec_reshape;
mov2 = mov_ori_reshape;
corr_pix=zeros(size(mov1,1),size(mov1,2));

for idim1=1:size(mov1,1)
    for idim2 = 1:size(mov1,2)
        seq1 = squeeze(mov1(idim1,idim2,tidx));
        seq2 = squeeze(mov2(idim1,idim2,tidx));
        corr_pix(idim1,idim2) = corr(seq1,seq2);
    end
end

figure;
imagesc((corr_pix));axis image;   ylim(xlim1);xlim(ylim1);%colormap gray

%% Analysis at different spatial frequencies ..
tidx=1:T;%400:800;
s1=40;s2=40;

mov_show = gather((Tkerinv*(Xdecode)')'); 
mov_show = temporal_filter(mov_show,ttf,0.7);
[mov_show1,xlim1,ylim1]= movie_reshape(mov_show,s1,s2,mask);
mov_show1 = mov_show1(xlim1(1)+2:xlim1(2)-2,ylim1(1)+2:ylim1(2)-2,:);


s1=40;s2=40;
mov_show = gather((Tkerinv*Xref')'); 
mov_show = temporal_filter(mov_show,ttf,0.7);
mov_show2= movie_reshape(mov_show,s1,s2,mask);
mov_show2 = mov_show2(xlim1(1)+2:xlim1(2)-2,ylim1(1)+2:ylim1(2)-2,:);

%mov_show2 = mov_show2*norm(mov_show1(:))/norm(mov_show2(:));

mask_sq = reshape(mask,[40,40])';

ss1=size(mov_show1,1);
ss2=size(mov_show1,2);
mov_err = mov_show1-mov_show2;
mov_err_fft_spatial = zeros(ss1,ss2);
for itime=tidx
mov_err_fft_spatial = mov_err_fft_spatial +abs(fft2(mov_err(:,:,itime)));
end
mov_err_fft_spatial =mov_err_fft_spatial /numel(tidx);

mov_show2_fft_spatial = zeros(ss1,ss2);
for itime=tidx
mov_show2_fft_spatial = mov_show2_fft_spatial +abs(fft2(mov_show2(:,:,itime)));
end
mov_show2_fft_spatial =mov_show2_fft_spatial /numel(tidx);


figure;
subplot(1,3,1);
imagesc(mov_err_fft_spatial);colormap gray;axis image
subplot(1,3,2);
imagesc(mov_show2_fft_spatial);colormap gray;axis image
subplot(1,3,3);
err_f_ratio = mov_err_fft_spatial./mov_show2_fft_spatial;
imagesc(err_f_ratio);colormap gray;axis image

sz1 = size(err_f_ratio,1);sz2 = size(err_f_ratio,2);
f_dist_mat = zeros(sz1,sz2);
for ix = 0:sz1-1
    for iy= 0:sz2-1
    ifx = 2*pi*ix/sz1;
    ify = 2*pi*iy/sz2;
    f_dist_mat(ix+1,iy+1) = norm([ifx,ify]); 
    end
end

figure;
a=f_dist_mat(:);
b=err_f_ratio(:);
[asort,idx] = sort(a,'ascend');
%plot(asort,b(idx));
%hold on;

a2=f_dist_mat(:);
b2=mov_show2_fft_spatial(:);
[a2sort,idx] = sort(a2,'ascend');
plotyy(asort,b(idx),a2sort,b2(idx));


