%% 
mov_ori = load('/Volumes/Analysis/nishal/original/movie.mat');
mov_ori_s2 = load('/Volumes/Analysis/nishal/original/mov_s2.mat');
mov_mod = load('/Volumes/Analysis/nishal/modified/movie.mat');
mov_s2 = load('/Volumes/Analysis/nishal/modified/mov_s2.mat');
mov_con_s10 = load('/Volumes/Analysis/nishal/modified/mov_con_s10.mat');

figure
for imovie=1:1200
   subplot(3,2,1);
   imagesc(mov_ori.mov(:,:,imovie));
   colormap(gray)
   caxis([0,255]);
   colorbar
   axis image
   title('Original')
   
   subplot(3,2,2);
   imagesc(mov_ori_s2.mov(:,:,imovie));
   colormap(gray)
   caxis([0,255]);
   colorbar
   axis image
   title('Original , s2')
   
    subplot(3,2,3);
   imagesc(mov_mod.mov(:,:,imovie));
   colormap(gray)
   caxis([0,255]);
   colorbar
   axis image
   title('Modified ')
   
    subplot(3,2,4);
   imagesc(mov_s2.mov(:,:,imovie));
   colormap(gray)
   caxis([0,255]);
   colorbar
   axis image
   title('Mod s2')
   
    subplot(3,2,5);
   imagesc(mov_con_s10.mov(:,:,imovie));
   colormap(gray)
   caxis([0,255]);
   colorbar
   axis image
   title('Mov con s10')
    
   pause(1/120);
end

