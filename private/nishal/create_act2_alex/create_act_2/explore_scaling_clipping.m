global_vars
toround=0;

if(toround==0)
cell_resp_orig=Ax(stas,mov_orignial,movie_time,n_cell);
cell_resp_null=Ax(stas,mov_modify_new,movie_time,n_cell);
else
cell_resp_orig=Ax(stas,mov_orignial,movie_time,n_cell);
cell_resp_null=Ax(stas,round(mov_modify_new),movie_time,n_cell);
end

%%
figure;
subplot(3,1,1);
hist(mov_orignial(:)+255/2,20);xlim([0,255]);
title('Origninal movie');

subplot(3,1,2);hist(mov_modify_new(:)+255/2,20);xlim([0,255])
title('modified, scaled, unclipped');
mm_new = mov_modify_new*1.3;
mm_new(mm_new>127.5)=127.5;
mm_new(mm_new<-127.5)=-127.5;
if(toround==1)

mm_new=round(mm_new);
end

subplot(3,1,3);hist(mm_new(:)+255/2,20);xlim([0,255])
cell_resp_scale=Ax(stas,mm_new,movie_time,n_cell);
title('modified,scaled,clipped');
%%

figure
subplot(5,1,1);
plot(cell_resp_orig);
ylim([-200,200]);
title('Original movie');



subplot(5,1,2);
plot(cell_resp_null);
ylim([-200,200])
title('Null movie response (bigger scale)');

subplot(5,1,3);
plot(cell_resp_null);
title('Null movie response (normal scale)');

subplot(5,1,4);
plot(cell_resp_scale);
ylim([-200,200]);
title('Modified movie, but scaled and clipped (bigger scale)');

subplot(5,1,5);
plot(cell_resp_scale);

title('Modified movie, but scaled and clipped (normal scale)');

%%

figure
for itime=1:100:movie_time
    subplot(2,2,1);
    imagesc(mov_orignial(:,:,itime));
    axis image
    colorbar
    colormap(gray)
   % title(sprintf('Time: %d',itime));
    caxis([-127.5,127.5]);
    
    subplot(2,2,2);
    imagesc(mov_modify_new(:,:,itime));
    axis image
    colorbar
   % title(sprintf('Time: %d',itime));
    colormap(gray)
    caxis([-127.5,127.5]);
        
    subplot(2,2,3);
    imagesc(mm_new(:,:,itime));
    axis image
    colorbar
  %  title(sprintf('Time: %d',itime));
    colormap(gray)
   caxis([-127.5,127.5]);
        
        
    pause(1/240)

end
