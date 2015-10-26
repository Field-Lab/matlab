
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


%% visualize the movie 
totalmask = reshape(sum(onpar.totalMaskAccept_log,2) + sum(offpar.totalMaskAccept_log,2),[32,16])';

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
    
    subplot(2,1,2);
    frame2 = reshape(original_tf(:,itime),[32,16])';
    frame2 = imfilter(frame2,imsm);
    imagesc(frame2.*totalmask);
    axis image;
    colormap gray;
    caxis([min(original_tf(:))*1,max(original_tf(:))*1]);
    title('Shown movie');
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

