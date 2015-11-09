
onpar = load('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2005_04_26_2/data002/ON_parasol.mat');
offpar = load('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2005_04_26_2/data002/OFF_parasol.mat');

%% calculate STAs

ncell_type1 = size(onpar.Y,1);
ncell_type2 = size(offpar.Y,1);

stas_celltype1 = cell(ncell_type1,1);
stas_celltype2 = cell(ncell_type2,1);

for icell=1:ncell_type1
stas_celltype1{icell} = onpar.maskedMovdd*onpar.Y(icell,:)'/sum(onpar.Y(icell,:));
stas_celltype1{icell} = stas_celltype1{icell}*onpar.ttf_log(:,icell)';
end


for icell=1:ncell_type2
stas_celltype2{icell} = offpar.maskedMovdd*offpar.Y(icell,:)'/sum(offpar.Y(icell,:));
stas_celltype2{icell} = stas_celltype2{icell}*offpar.ttf_log(:,icell)';
end

%% reconstruct 

mov_recons1 = zeros(size(onpar.maskedMovdd,1),size(onpar.maskedMovdd,2));

idx = 1:size(onpar.maskedMovdd,2);

for icell=1:ncell_type1
    icell
    for ispk=idx(onpar.Y(icell,:)>0 & idx>30);
    mov_recons1(:,ispk-29:ispk) = mov_recons1(:,ispk-29:ispk) + stas_celltype1{icell}(:,end:-1:1);
    end    
end

mov=mov_recons1;
mov_recons1_filt= filterMov_cone(mov,ones(size(mov,1),1),squeeze(onpar.ttf_avg));


mov_recons2 = zeros(size(onpar.maskedMovdd,1),size(onpar.maskedMovdd,2));

idx = 1:size(offpar.maskedMovdd,2);

for icell=1:ncell_type2
    icell
    for ispk=idx(offpar.Y(icell,:)>0 & idx>30);
    mov_recons2(:,ispk-29:ispk) = mov_recons2(:,ispk-29:ispk) + stas_celltype2{icell}(:,end:-1:1);
    end    
end

mov=mov_recons2;
mov_recons2_filt= filterMov_cone(mov,ones(size(mov,1),1),squeeze(offpar.ttf_avg));

mov_recons = mov_recons1+mov_recons2;
mov=mov_recons;
mov_recons_filt_on  = filterMov_cone(mov,ones(size(mov,1),1),squeeze(onpar.ttf_avg));
mov_recons_filt_off  = filterMov_cone(mov,ones(size(mov,1),1),squeeze(offpar.ttf_avg));


%% visualize the movie - flat movies

totalmask = reshape(sum(onpar.totalMaskAccept_log,2) + sum(offpar.totalMaskAccept_log,2),[32,16])';

[r,c] = find(totalmask>0);
r_range = [min(r)+3,max(r)-3];
c_range = [min(c)+3,max(c)-3];
imsm = ones(4,4)/16;
mov = gather(mov_recons_filt_off);
original_tf = gather(offpar.maskedMovdd);

figure;
for itime =100:1000; 
    subplot(2,1,1);
    frame1 = reshape(mov(:,itime),[32,16])';
   % frame1 = imfilter(frame1,imsm);
    imagesc(frame1.*totalmask);
    axis image;
    colormap gray;
    caxis([min(mov(:))*1,max(mov(:))*1]);
    title('Reconstructed');
    xlim(c_range);ylim(r_range);
    
    subplot(2,1,2);
    frame2 = reshape(original_tf(:,itime),[32,16])';
    frame2 = imfilter(frame2,imsm);
    imagesc(frame2.*totalmask);
    axis image;
    colormap gray;
    caxis([min(original_tf(:))*1,max(original_tf(:))*1]);
    title('Shown movie');
       xlim(c_range);ylim(r_range);
    pause%(1/50);
    %title(sprintf('%d',itime));
end


%% fft mask

fft_abs = zeros(16,32); 
for iframe = 1:size(mov,3)
    fr = reshape(mov(:,itime),[32,16])';
    fft_abs = fft_abs+ abs(fft2(fr));
end

fft_abs = fft_abs/size(mov,2);
figure;
contour(log(fft_abs),10)

%%  temporal frequency filtering

temp_f_cutoff=0.7;

T = size(onpar.maskedMovdd,2);
figure('Color','w');
ftt_off = fft(offpar.ttf_avg,T);
plot([0:2/T:2*(T-1)/T],abs(ftt_off))
hold on;

ftt_on = fft(onpar.ttf_avg,T);
plot([0:2/T:2*(T-1)/T],abs(ftt_on))
hold on;

ftt_total = abs(ftt_on)+abs(ftt_off);
plot([0:2/T:2*(T-1)/T],ftt_total/2);
hold on;

fft_mask = double(ftt_total > temp_f_cutoff*max(ftt_total));
plot([0:2/T:2*(T-1)/T],fft_mask);
legend('Off Parasol','On Parasol','Average','Selected mask (0.9>max)');


% filter reconstructed movie
movin = mov_recons;
movout = 0*movin;
for idim = 1:size(movin,1)
   movout(idim,:)= ifft(fft(movin(idim,:),T).*fft_mask');
end

movout = movout * norm(movin(:)) / norm(movout(:));
mov_recons_fft = movout;


% filter original movie
movin = onpar.mov;movout= 0*movin;
for idim1 = 1:size(movin,1)
    for idim2 = 1:size(movin,2)
    movout(idim1,idim2,:) = ifft(fft(squeeze(movin(idim1,idim2,:)),T).*fft_mask);
    end
end
movout = movout * norm(movin(:)) / norm(movout(:));
movout = movout *norm(mov_recons_fft(:))/norm(movout(:));
mov_orig_fft = movout;


%% spatial filtering
spatial_f_cutoff = 0.5;

% find frequencies relevant for STA
fft_xx = zeros(32,16);
stas =stas_celltype1;
for icell=1:length(stas)
xx = reshape(stas{icell}(:,4),[32,16]);
fft_xx = fft_xx + abs(fft2(xx));
end

stas =stas_celltype2;
for icell=1:length(stas)
xx = reshape(stas{icell}(:,4),[32,16]);
fft_xx = fft_xx + abs(fft2(xx));
end

fft_sp_mask = double(fft_xx>=spatial_f_cutoff*max(fft_xx(:)));
figure('Color','w');
subplot(2,1,1);
imagesc(fft_xx');
axis image
title('Average Fourier tranform(magnitude) of RFs')
subplot(2,1,2);
imagesc(fft_sp_mask')
axis image 
title('Selected frequencies')

 movin = mov_recons_fft;movout =zeros(32,16,T);
 for itime = 1:T
 xx =  reshape(movin(:,itime),[32,16]);
 movout(:,:,itime) = ifft2(fft2(xx).*fft_sp_mask);
 end
 movout = movout * norm(movin(:)) / norm(movout(:));
 mov_recons_fft2= movout;
 
  movin = mov_orig_fft;movout =zeros(32,16,T);
 for itime = 1:T
 xx =  movin(:,:,itime);
 movout(:,:,itime) = ifft2(fft2(xx).*fft_sp_mask);
 end
 movout = movout * norm(movin(:)) / norm(movout(:));
 mov_orig_fft2= movout;
%% visualize the movie - flat movies

totalmask = reshape(sum(onpar.totalMaskAccept_log,2) + sum(offpar.totalMaskAccept_log,2),[32,16])';

[r,c] = find(totalmask>0);
r_range = [min(r)+3,max(r)-3];
c_range = [min(c)+3,max(c)-3];
    
mov =mov_recons_fft2;
original_tf = mov_orig_fft2;

figure;
for itime =206:1000; 
    subplot(2,1,1);
   frame1= mov(:,:,itime)'; 
   imagesc(frame1.*totalmask);
    axis image;
    colormap gray;
    %caxis([min(mov(:))*1,max(mov(:))*1]);
    title('Reconstructed');
    xlim(c_range);ylim(r_range);
    set(gca,'xTick',[]);
    set(gca,'yTick',[]);
    
    subplot(2,1,2);
    frame2 = original_tf(:,:,itime)';%reshape(original_tf(:,itime),[32,16])';
    imagesc(frame2.*totalmask);
    axis image;
    colormap gray;
    %caxis([min(original_tf(:))*1,max(original_tf(:))*1]);
    title('Shown movie');
    set(gca,'xTick',[]);
    set(gca,'yTick',[]);
    
       xlim(c_range);ylim(r_range);
    suptitle(sprintf('%d',itime));
       pause
end

%% Find correlation between pixel time courses

mov1 = mov_recons_fft2;
mov2 = mov_orig_fft2;
corr_pix=zeros(size(mov1,1),size(mov1,2));

for idim1=1:size(mov1,1)
    for idim2 = 1:size(mov1,2)
        seq1 = squeeze(mov1(idim1,idim2,:));
        seq2 = squeeze(mov2(idim1,idim2,:));
        corr_pix(idim1,idim2) = corr(seq1,seq2);
    end
end

figure;
imagesc(corr_pix');
axis image
%%
mov_recons_fft2_short = mov_recons_fft2(:,:,1:2000);
mov_orig_fft2_short = mov_orig_fft2(:,:,1000:2000);
%% show time courses 
mov1 = mov_recons_fft2_short;
mov2 = mov_orig_fft2_short;

time_idx = [120:241];
pix_list =[18,11;
            23,8;
            9,8;
            10,12;
            14,6];

jump=0.4;

h=figure('Color','w');
for ipix = 1:size(pix_list,1) 
    idim1=pix_list(ipix,1) %1:size(mov1,1)
     idim2 = pix_list(ipix,2) %1:size(mov1,2)
        seq1 = squeeze(mov1(idim1,idim2,time_idx)); seq1=seq1/norm(seq1);
        seq2 = squeeze(mov2(idim1,idim2,time_idx)); seq2 =seq2/norm(seq2);
        
        plot([0:length(time_idx)-1]/120,seq1-(ipix-1)*jump,'r');
        hold on;
        plot([0:length(time_idx)-1]/120,seq2-(ipix-1)*jump,'k')
        hold on;
        set(gca,'yTick',[]);
        hold on;
        text(max((length(time_idx)-1)/120),-(ipix-1)*jump,sprintf('%0.03f',corr_pix(idim1,idim2)));
       % title(sprintf('%0.02f',corr_pix(idim1,idim2)));
end

legend('Reconstruction','Input stimulus');


%% Reconstruction error frequency analysis .. 



