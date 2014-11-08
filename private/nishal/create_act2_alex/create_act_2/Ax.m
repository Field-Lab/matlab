function cell_resp = Ax(stas,mov,movie_time,n_cell)
% Nishal P. shah
% August 2014


%KKt=filt_mat*filt_mat';
filt_dim1=size(mov,1);
filt_dim2=size(mov,2);
filt_len=size(stas{1},4);

mov2=zeros(filt_dim1,filt_dim2,movie_time+filt_len-1);
mov2(:,:,filt_len:movie_time+filt_len-1)=mov; % Append zeros before the movie

cell_resp=zeros(movie_time,n_cell);
sz=max(size(mov2,3)-size(stas{1},4) + 1, 0);
parfor icell=1:n_cell
cell_resp(:,icell)=reshape(convn(mov2,stas{icell}(end:-1:1,end:-1:1,:),'valid'),[sz,1]);
end

end