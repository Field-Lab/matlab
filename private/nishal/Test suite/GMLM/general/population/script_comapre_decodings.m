%% Here , we compare two decodings

%% laod decodings 1
% 
% decode1 = load('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_ifit1.mat');
% decode17 = load('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_ifit17.mat');
% 
% Xdecode_a = decode1.Xdecode;
% Xdecode_b = decode17.Xdecode;

%% load decodings 2
dec = load('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_ifit_1_17_T4000.mat');
Xdecode_a = dec.Xdecode_ifit1;
Xdecode_b = dec.Xdecode_ifit17;


%%
tidx=1:T;%400:800;
s1=40;s2=40;

spatial_cutoff_u=1.2;
spatial_cutoff_l=0.25;


mov_show = gather((Tkerinv*(Xdecode_a)')'); 
mov_show = temporal_filter(mov_show,ttf,0.7);
[mov_show,xlim1,ylim1]= movie_reshape(mov_show,s1,s2,mask);
mov_dec_reshape_a = spatial_filter(mov_show,stas_celltype1,spatial_cutoff_l,spatial_cutoff_u,s1,s2);
mov_dec_cut_a = mov_dec_reshape_a(xlim1(1):xlim1(2),ylim1(1):ylim1(2),:);


mov_show = gather((Tkerinv*(Xdecode_b)')'); 
mov_show = temporal_filter(mov_show,ttf,0.7);
[mov_show,xlim1,ylim1]= movie_reshape(mov_show,s1,s2,mask);
mov_dec_reshape_b = spatial_filter(mov_show,stas_celltype1,spatial_cutoff_l,spatial_cutoff_u,s1,s2);
mov_dec_cut_b = mov_dec_reshape_b(xlim1(1):xlim1(2),ylim1(1):ylim1(2),:);



s1=40;s2=40;
mov_show = gather((Tkerinv*Xref')'); 
mov_show = temporal_filter(mov_show,ttf,0.7);
mov_show= movie_reshape(mov_show,s1,s2,mask);
mov_ori_reshape = spatial_filter(mov_show,stas_celltype1,spatial_cutoff_l,spatial_cutoff_u,s1,s2);
mov_ori_cut = mov_ori_reshape(xlim1(1):xlim1(2),ylim1(1):ylim1(2),:);

figure;
icnt=1;
for itime=[1069,1170,1198]%size(Xref,2)
   
    subplot(3,3,icnt);
   imagesc(repelem(mov_dec_cut_b(:,:,itime),20,20));axis image;colormap gray;%caxis([min(mov_dec_reshape(:)), 0.5*max(mov_dec_reshape(:))]);
   set(gca,'xTick',[]);set(gca,'yTick',[]);
   icnt=icnt+1; 
   
   subplot(3,3,icnt);
   imagesc(repelem(mov_dec_cut_a(:,:,itime),20,20));axis image;colormap gray;%caxis([min(mov_dec_reshape(:)), 0.5*max(mov_dec_reshape(:))]);
   set(gca,'xTick',[]);set(gca,'yTick',[]);
   icnt=icnt+1;
   
   subplot(3,3,icnt);
   imagesc(repelem(mov_ori_cut(:,:,itime),20,20));axis image;colormap gray;%caxis([min(mov_ori_reshape(:)),0.5*max(mov_ori_reshape(:))]);
    set(gca,'xTick',[]);set(gca,'yTick',[]);
   icnt=icnt+1; 
   
   itime
  
end

%% plot time courses
figure;

subplot(3,1,1);
tidx=2000:2400;
ipix = 3; jpix=8;
plot(tidx/120,squeeze(mov_ori_cut(ipix,jpix,tidx)));
hold on;
scale=norm(squeeze(mov_ori_cut(ipix,jpix,tidx)))/norm(squeeze(mov_dec_cut_a(ipix,jpix,tidx)));
plot(tidx/120,squeeze(scale*mov_dec_cut_a(ipix,jpix,tidx)),'LineWidth',2);
hold on;
scale=norm(squeeze(mov_ori_cut(ipix,jpix,tidx)))/norm(squeeze(mov_dec_cut_b(ipix,jpix,tidx)));
plot(tidx/120,squeeze(scale*mov_dec_cut_b(ipix,jpix,tidx)));
xlim([tidx(1),tidx(end)]/120);
cr =corr(squeeze(mov_ori_cut(ipix,jpix,1000:3000)),squeeze(mov_dec_cut_a(ipix,jpix,1000:3000)))
title(sprintf('correlation %02f',cr));
set(gca,'yTick',[]);

subplot(3,1,2);
tidx=1000:1400;
ipix =5; jpix=4;
plot(tidx/120,squeeze(mov_ori_cut(ipix,jpix,tidx)));
hold on;
scale=norm(squeeze(mov_ori_cut(ipix,jpix,tidx)))/norm(squeeze(mov_dec_cut_a(ipix,jpix,tidx)));
plot(tidx/120,squeeze(scale*mov_dec_cut_a(ipix,jpix,tidx)),'LineWidth',2);
hold on;
scale=norm(squeeze(mov_ori_cut(ipix,jpix,tidx)))/norm(squeeze(mov_dec_cut_b(ipix,jpix,tidx)));
plot(tidx/120,squeeze(scale*mov_dec_cut_b(ipix,jpix,tidx)));
xlim([tidx(1),tidx(end)]/120);
cr =corr(squeeze(mov_ori_cut(ipix,jpix,1000:3000)),squeeze(mov_dec_cut_a(ipix,jpix,1000:3000)))
title(sprintf('correlation %02f',cr));
set(gca,'yTick',[]);


subplot(3,1,3);
tidx=2500:2900;
ipix = 6; jpix=9;
plot(tidx/120,squeeze(mov_ori_cut(ipix,jpix,tidx)));
hold on;
scale=norm(squeeze(mov_ori_cut(ipix,jpix,tidx)))/norm(squeeze(mov_dec_cut_a(ipix,jpix,tidx)));
plot(tidx/120,squeeze(scale*mov_dec_cut_a(ipix,jpix,tidx)),'LineWidth',2);
hold on;
scale=norm(squeeze(mov_ori_cut(ipix,jpix,tidx)))/norm(squeeze(mov_dec_cut_b(ipix,jpix,tidx)));
plot(tidx/120,squeeze(scale*mov_dec_cut_b(ipix,jpix,tidx)));
xlim([tidx(1),tidx(end)]/120);
corr(squeeze(mov_ori_cut(ipix,jpix,1000:3000)),squeeze(mov_dec_cut_a(ipix,jpix,1000:3000)))
cr =corr(squeeze(mov_ori_cut(ipix,jpix,1000:3000)),squeeze(mov_dec_cut_a(ipix,jpix,1000:3000)))
set(gca,'yTick',[]);
title(sprintf('correlation %02f',cr));
legend('True stimulus','Decoded by 1 SU','Decoded by 17 SU');

xlabel('Time (s)')

%% Outline of STAs

tm = total_mask_log(:,cellsChoose);

ncell = size(tm,2);
for icell=1:ncell
xx = reshape(tm(:,icell),[s1,s2])';
xx=xx(xlim1(1):xlim1(2),ylim1(1):ylim1(2));
ttm{icell} = xx;
end

B_use = cell(ncell,1);
for icell=1:ncell
 B = bwboundaries(repelem(ttm{icell},20,20));
B_use{icell} = B{1};
end
cols = distinguishable_colors(ncell+12); cols = cols(5:end,:);

    
%% Correlations
mask2D = reshape(mask,[40,40])';
mask2D = mask2D(xlim1(1):xlim1(2),ylim1(1):ylim1(2));

% dec a and ori
mov1 = mov_dec_cut_a;
mov2 = mov_ori_cut;
tidx = 1000:3000;

corr_a_ori  = corr_mov(mov1(:,:,tidx),mov2(:,:,tidx));
h=figure;imagesc(repelem(corr_a_ori.*mask2D,20,20));axis image;colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]);  % ylim(xlim1);xlim(ylim1);%colormap gray
title('a and ori');
plot_cell_bdry(B_use,cols,1,1)
%print(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_figures/corr_a_ori.pdf','-dpdf');

% dec b and ori
mov1 = mov_dec_cut_b;
mov2 = mov_ori_cut;

corr_b_ori  = corr_mov(mov1(:,:,tidx),mov2(:,:,tidx));
h=figure;imagesc(repelem(corr_b_ori.*mask2D,20,20));axis image;colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]); % ylim(xlim1);xlim(ylim1);%colormap gray
title('b and ori');
plot_cell_bdry(B_use,cols,1,1)
%print(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_figures/corr_b_ori.pdf','-dpdf');


corr_ba_ori_difff = corr_b_ori-corr_a_ori;
h=figure; imagesc(repelem(corr_ba_ori_difff.*mask2D,20,20)); axis image;colormap gray%caxis([-max(abs(corr_ba_ori_diff(:)))*0.5, max(abs(corr_ba_ori_diff(:)))*0.5])
colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]);
title('Difference of corr between two movies');
plot_cell_bdry(B_use,cols,1,1)
%print(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_figures/corr_b_ori_minus_corr_a_ori.pdf','-dpdf');


% maximal difference between the two reconstructions
mov_diff = mov_dec_cut_a - mov_dec_cut_b;
diff_mag = squeeze(sqrt(sum(sum(mov_diff.^2,1),2)));
diff_mag = diff_mag.*double([1:T]>=400 & [1:T]<=3600)';
t_chosen = diff_mag>prctile(diff_mag(tidx),90);

corr_a_ori_diff  = corr_mov(mov_dec_cut_a(:,:,t_chosen),mov_ori_cut(:,:,t_chosen));
h=figure;imagesc(repelem(corr_a_ori_diff.*mask2D,20,20));axis image; colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]); % ylim(xlim1);xlim(ylim1);%colormap gray
title('corr a and ori at max diff with b');
plot_cell_bdry(B_use,cols,1,1)
%print(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_figures/corr_a_ori_diff_withb.pdf','-dpdf');


corr_b_ori_diff  = corr_mov(mov_dec_cut_b(:,:,t_chosen),mov_ori_cut(:,:,t_chosen));
h=figure;imagesc(repelem(corr_b_ori_diff.*mask2D,20,20));axis image; colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]);   % ylim(xlim1);xlim(ylim1);%colormap gray
title('corr b and ori at max diff with a');
plot_cell_bdry(B_use,cols,1,1)
%print(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_figures/corr_b_ori_diff_witha.pdf','-dpdf');


corr_ba_ori_diff = corr_b_ori_diff-corr_a_ori_diff;
h=figure; imagesc(repelem(corr_ba_ori_diff.*mask2D,20,20)); axis image;colormap gray;%caxis([-max(abs(corr_ba_ori_diff(:)))*0.5, max(abs(corr_ba_ori_diff(:)))*0.5])
colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]);
title('Difference of corr between two movies');
plot_cell_bdry(B_use,cols,1,1)
%print(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_figures/corr_b_ori_minus_corr_a_ori__maximal_diff_b_a.pdf','-dpdf');


h=figure;plot(corr_a_ori_diff(:),corr_b_ori_diff(:),'.');hold on; plot([0,1],[0,1],'g');
title('scatter plot between pixel correlations for a and b');


% pixelwise maximal difference!
pc=90;
[corr_pixa,corr_pixb] = corr_mov_maximal_diff_pixelwise(mov_dec_cut_a(:,:,tidx),mov_dec_cut_b(:,:,tidx),mov_ori_cut(:,:,tidx),pc);

h=figure;imagesc(repelem(corr_pixa.*mask2D,20,20));axis image;  % ylim(xlim1);xlim(ylim1);%colormap gray
colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]);
title('a and ori pixelwise maximal difference');
plot_cell_bdry(B_use,cols,1,1)
%print(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_figures/corr_a_ori_pixelwise_maximal_diff_b.pdf','-dpdf');


h=figure;imagesc(repelem(corr_pixb.*mask2D,20,20));axis image;  % ylim(xlim1);xlim(ylim1);%colormap gray
colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]);
title('b and ori pixelwise maximal difference');
plot_cell_bdry(B_use,cols,1,1)
%print(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_figures/corr_b_ori_pixelwise_maximal_diff_a.pdf','-dpdf');

h=figure;imagesc(repelem((corr_pixb-corr_pixa).*mask2D,20,20));axis image;  % ylim(xlim1);xlim(ylim1);%colormap gray
title('b-a and ori pixelwise maximal difference');colormap gray
colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]);
plot_cell_bdry(B_use,cols,1,1)
%print(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_figures/corr_b_ori_minus_corr_a_ori_pixelwise_maximal_diff_ab.pdf','-dpdf');


h=figure;plot(corr_pixa(:),corr_pixb(:),'.');hold on; plot([0,1],[0,1],'g');colormap gray
title('scatter plot between pixel correlations for a and b');


% when decoding differs from original, where is it the worst? 
mov_diff = mov_dec_cut_a - mov_ori_cut;
diff_mag = squeeze(sqrt(sum(sum(mov_diff.^2,1),2)));
diff_mag = diff_mag.*double([1:T]>=1000 & [1:T]<=3000)';
t_chosen = diff_mag>prctile(diff_mag(tidx),90);

corr_a_ori2  = corr_mov(mov_dec_cut_a(:,:,t_chosen),mov_ori_cut(:,:,t_chosen));
h=figure;imagesc(repelem(corr_a_ori2.*mask2D,20,20));axis image;  % ylim(xlim1);xlim(ylim1);%colormap gray
colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]);
title('a and ori differ max overall');
plot_cell_bdry(B_use,cols,1,1)
%print(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_figures/corr_a_ori_diff_ORI.pdf','-dpdf');


mov_diff = mov_dec_cut_b - mov_ori_cut;
diff_mag = squeeze(sqrt(sum(sum(mov_diff.^2,1),2)));
diff_mag = diff_mag.*double([1:T]>=1000 & [1:T]<=3000)';
t_chosen = diff_mag>prctile(diff_mag(tidx),90);

corr_b_ori2 = corr_mov(mov_dec_cut_b(:,:,t_chosen),mov_ori_cut(:,:,t_chosen));
h=figure;imagesc(repelem(corr_b_ori2.*mask2D,20,20));axis image;  % ylim(xlim1);xlim(ylim1);%colormap gray
colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]);
title('b and ori differ max overall');
plot_cell_bdry(B_use,cols,1,1)
%print(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_figures/corr_b_ori_diff_ORI.pdf','-dpdf');


h=figure;imagesc(repelem((corr_b_ori2-corr_a_ori2).*mask2D,20,20));axis image;  % ylim(xlim1);xlim(ylim1);%colormap gray
colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]);
title('b-a and x and ori differ max overall');colormap gray
plot_cell_bdry(B_use,cols,1,1)
%print(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/decode_figures/corr_b_ori_minus_corr_a_ori_diff_ORI.pdf','-dpdf');

