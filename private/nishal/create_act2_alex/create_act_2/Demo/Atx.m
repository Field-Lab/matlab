function mov_new=Atx(stas,cell_resp,filt_dim1,filt_dim2,filt_len,movie_time,n_cell)

mov_new=zeros(filt_dim1,filt_dim2,filt_len+movie_time-1);
parfor icell=1:n_cell
mov_new=mov_new+convn(stas{icell}(:,:,end:-1:1),reshape(cell_resp(:,icell),[1,1,length(cell_resp(:,icell))])); % DEBUG?
end
mov_new=mov_new(:,:,filt_len:filt_len+movie_time-1);

end