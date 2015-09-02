function display_movies(recon_movie, original_movie)

pause_length = 0.01;
movie_length = length(original_movie);

figure;

for i_frame = 1:movie_length
    
subplot(2,1,1)
imagesc(original_movie(:,:,i_frame)')
colormap gray
caxis([0 1])
axis image
axis off

subplot(2,1,2)
imagesc(recon_movie(:,:,i_frame)')
colormap gray
caxis([0 1])
axis image
axis off

pause(pause_length)

end



end