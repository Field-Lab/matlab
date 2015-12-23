tic;
filt_len=size(stas{1},4);

% making sparse toeplitz ?filt_mat= sparse([],[],[],movie_time*n_cell,filt_dim1*filt_dim2*movie_time,n_cell*movie_time*6*320*160);
%filt_mat= sparse([],[],[],movie_time*1,filt_dim1*filt_dim2*movie_time,n_cell*movie_time*filt_len*filt_dim1*filt_dim2);
%filt_mat=cell(n_cell,1);
%inv_filt_mat=cell(n_cell,1);

nzero=n_cell*movie_time*filt_len*filt_dim1*filt_dim2;

x_log=zeros(nzero,1);
y_log=zeros(nzero,1);
v_log=zeros(nzero,1);
cnt=0;

for icell=1:n_cell
k=stas{icell};


icell
for itime=1:movie_time
   
    for iitime=1:min(itime,filt_len)
    k_t=k(:,:,iitime);   
    k_t=k_t(:)';
    %filt_mat(itime+(icell-1)*(movie_time),itime*filt_dim1*filt_dim2-(iitime)*filt_dim1*filt_dim2+1:itime*filt_dim1*filt_dim2-(iitime-1)*filt_dim1*filt_dim2)=k_t(:)';
    for ispace=1:filt_dim1*filt_dim2
        cnt=cnt+1;
    x_log(cnt)=itime+(icell-1)*(movie_time);
    y_log(cnt)=itime*filt_dim1*filt_dim2-(iitime)*filt_dim1*filt_dim2+ispace;
    v_log(cnt)=k_t(ispace);
        
    end
    
    end
end


end
x_log=x_log(1:cnt);
y_log=y_log(1:cnt);
v_log=v_log(1:cnt);

filt_mat=sparse(x_log,y_log,v_log,movie_time*n_cell,filt_dim1*filt_dim2*movie_time,n_cell*movie_time*filt_len*filt_dim1*filt_dim2);
toc;
tic;
%inv_filt_mat=(filt_mat{icell}*filt_mat{icell}')\filt_mat{icell};
toc;
