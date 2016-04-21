spikes = load('/Volumes/Lab/Users/crhoades/Colleen/matlab/private/colleen/New Cell Types/Jitter/data026_cf_all_smooth.mat')
cellID =1;
stas = load('/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data026cf_spikes/Cell 481.mat')

figure;
for istep=15
    istep
imagesc(mean(stas.temp(:,:,:,istep),3)');axis image;colormap gray;
caxis([min(stas.temp(:)),max(stas.temp(:))]);
pause
end

xx = 200:350;yy = 100:200;
mask=zeros(640,320);
mask(xx,yy)=1;
mask=mask';

%% get movie

mm = zeros(sum(mask(:)),200*1250);
icnt=1;
for itime=1:1250
    mov = load(sprintf('/Volumes/Lab/Users/crhoades/JitterMovie/2016-02-17-6/data026/movie_block_%d.mat',itime));
    mmov = squeeze(mean(mov.current_movie,3));
    for iicnt=1:200
        xx=mmov(:,:,iicnt);
        mm(:,icnt)=xx(logical(mask));
        icnt=icnt+1;
    end
end

