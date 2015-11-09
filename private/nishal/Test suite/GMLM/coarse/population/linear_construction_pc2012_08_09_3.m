

onpar = load('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/On_par2.mat');
offpar = load('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/Off_par2.mat');

%% calculate STAs

ncell_type1 = size(onpar.binnedSpikeResponses_coll,1);
ncell_type2 = size(offpar.binnedSpikeResponses_coll,1);

stas_celltype1 = cell(ncell_type1,1);
stas_celltype2 = cell(ncell_type2,1);

for icell=1:ncell_type1
stas_celltype1{icell} = onpar.maskedMovdd*onpar.binnedSpikeResponses_coll(icell,:)'/sum(onpar.binnedSpikeResponses_coll(icell,:));
stas_celltype1{icell} = stas_celltype1{icell}*onpar.ttf';
end


for icell=1:ncell_type2
stas_celltype2{icell} = offpar.maskedMovdd*offpar.binnedSpikeResponses_coll(icell,:)'/sum(offpar.binnedSpikeResponses_coll(icell,:));
stas_celltype2{icell} = stas_celltype2{icell}*offpar.ttf';
end

%% reconstruct 

mov_recons = zeros(size(onpar.maskedMovdd,1),size(onpar.maskedMovdd,2));

idx = 1:size(onpar.maskedMovdd,2);

for icell=1:ncell_type1
    icell
    for ispk=idx(onpar.binnedSpikeResponses_coll(icell,:)>0 & idx>30);
    mov_recons(:,ispk-29:ispk) = mov_recons(:,ispk-29:ispk) + stas_celltype1{icell}(:,end:-1:1);
    end    
end

idx = 1:size(offpar.maskedMovdd,2);

for icell=1:ncell_type2
    icell
    for ispk=idx(offpar.binnedSpikeResponses_coll(icell,:)>0 & idx>30);
    mov_recons(:,ispk-29:ispk) = mov_recons(:,ispk-29:ispk) + stas_celltype2{icell}(:,end:-1:1);
    end    
end

figure;
for itime =50:1000; 
    subplot(2,1,1);
    imagesc(reshape(mov_recons(:,itime),[80,40])');
    axis image;
    colormap gray;
    
    subplot(2,1,2);
    imagesc(reshape(onpar.maskedMovdd(:,itime),[80,40])');
    axis image;
    colormap gray;
    
    pause(1/50);
    title(sprintf('%d',itime));
end

