function [cMaps,h]=makecMap(destination_mat,movIDs,rawMovFrames)

noMovs=length(movIDs);

cMaps=cell(noMovs,1);
for imov=1:noMovs
[stim,height,width,header_size] = get_raw_movie(sprintf('%s/%d.rawMovie',destination_mat,movIDs(imov)),rawMovFrames,1);

subtract_movies=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies,[rawMovFrames,1,1]);
movie=movie/255;
condMovies=permute(movie,[2,3,1]);

cMaps{imov}  = sqrt(sum(condMovies.^2,3));
end

h=figure('Color','w');
for imov=1:noMovs
subplot(1,noMovs,imov);
imagesc(cMaps{imov});
colorbar
colormap gray
axis image
caxis([min(cMaps{2}(:)),max(cMaps{2}(:))+0.5]);

end
end