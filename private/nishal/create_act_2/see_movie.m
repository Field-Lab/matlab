imov=1;
%load(sprintf('/Volumes/Analysis/nishal/movies/mov_%d',imov));
figure;

for itime=1:100:movie_time
    subplot(2,2,1);
    imagesc(mov_log{imov}(:,:,itime));
    axis image
    colorbar
    colormap(gray)
    title(sprintf('Time: %d',itime));
    
    subplot(2,2,2);
    imagesc(mov_modify_new(:,:,itime));
    axis image
    colorbar
    title(sprintf('Time: %d',itime));
    colormap(gray)
        
    subplot(2,2,3);
    imagesc(mov_log{imov}(:,:,itime)-mov_modify_new(:,:,itime));
    axis image
    colorbar
    title(sprintf('Time: %d',itime));
    colormap(gray)
    
    pause(1/240)

end

%%
cell_resp_orig=Ax(stas,mov_log{imov},movie_time,n_cell);
cell_resp_null=Ax(stas,mov_modify_new,movie_time,n_cell);


figure
subplot(3,1,1);
plot(cell_resp_orig);
ylim([-200,200]);
title('Original movie');



subplot(3,1,2);
plot(cell_resp_null);
ylim([-200,200])
title('Null movie response (bigger scale)');

subplot(3,1,3);
plot(cell_resp_null);
title('Null movie response (normal scale)');

%%
clear var
figure;
subplot(2,1,1);
hist(var(cell_resp_orig),10)
title('Variance of Origional response');

subplot(2,1,2);
hist(var(cell_resp_null),10)
title('Variance of null response ');

%% pixel histogram

clear hist
figure;
subplot(2,1,1);
hist(mov_log{imov}(:),20);
title('Pixel histogram: Original movie');

subplot(2,1,2);
hist(mov_modify_new(:),20);
title('Pixel histogram: Modified movie');
