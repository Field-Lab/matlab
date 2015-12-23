q=zeros(size(stas{1},1),size(stas{1},2),size(stas{1},4));
for itime=1:30
    q(:,:,itime)=stas{1}(end:-1:1,end:-1:1,itime);
end

%%
% K*c
mov2=zeros(filt_dim1,filt_dim2,movie_time+filt_len-1);
mov2(:,:,filt_len:movie_time+filt_len-1)=mov;

q=stas{2}(end:-1:1,end:-1:1,:);

ax=convn(mov2,q,'valid');
ax=reshape(ax,[size(ax,3),1]);

%%
mov_reshape2=zeros(filt_dim1*filt_dim2*movie_time,1);
itime=0;
for i=1:filt_dim1*filt_dim2:length(mov_reshape2)
    itime=itime+1;
    xx=ax(:,:,itime);
  mov_reshape2(i:i+filt_dim1*filt_dim2-1)=  xx(:);
end

%%
qt = stas{2}(end:-1:1,end:-1:1,end:-1:1);
axx=reshape(ax,[1,1,length(ax)]);
mov_new=convn(qt,axx);

%%
%% 
% Do for all cells
mov2=zeros(filt_dim1,filt_dim2,movie_time+filt_len-1);
mov2(:,:,filt_len:movie_time+filt_len-1)=mov; % Append zeros before the movie

cell_resp=zeros(movie_time,n_cell);

for icell=1:n_cell
q=stas{icell}(end:-1:1,end:-1:1,:);
ax=convn(mov2,q,'valid');
ax=reshape(ax,[size(ax,3),1]);
cell_resp(:,icell)=ax;
end

cell_resp = Ax(stas,mov2,movie_time,n_cell);

mov_new=zeros(32,64,149);
for icell=1:n_cell
qt = stas{icell}(end:-1:1,end:-1:1,end:-1:1);
ax=cell_resp(:,icell);
axx=reshape(ax,[1,1,length(ax)]);
mov_new=mov_new+convn(qt,axx);
end
mov_new=mov_new(:,:,filt_len:filt_len+movie_time-1);

%%
figure;
for itime=1:120%size(mov_new,3)
imagesc(mov_new(:,:,itime));
colormap(gray);
title(sprintf('Time: %d:',itime));
pause(1/120)
end
%%
figure;
for itime=1:30%size(mov_new,3)
imagesc(qt(:,:,itime));
colormap(gray);
title(sprintf('Time: %d:',itime));
pause(0.5)
end




