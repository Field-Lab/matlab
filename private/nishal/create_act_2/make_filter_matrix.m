function [filt_mat,inv_filt_mat]=make_filter_matrix(k,movie_time,filt_len,filt_dim1,filt_dim2)

%k=stas{icell};
nzero=1*movie_time*filt_len*filt_dim1*filt_dim2;
x_log=zeros(nzero,1);
y_log=zeros(nzero,1);
v_log=zeros(nzero,1);
cnt=0;

for itime=1:movie_time
   
    for iitime=1:min(itime,filt_len)
    k_t=k(:,:,iitime);   
    k_t=k_t(:)';
    %filt_mat(itime+(icell-1)*(movie_time),itime*filt_dim1*filt_dim2-(iitime)*filt_dim1*filt_dim2+1:itime*filt_dim1*filt_dim2-(iitime-1)*filt_dim1*filt_dim2)=k_t(:)';
    for ispace=1:filt_dim1*filt_dim2
        cnt=cnt+1;
    x_log(cnt)=itime+(1-1)*(movie_time);
    y_log(cnt)=itime*filt_dim1*filt_dim2-(iitime)*filt_dim1*filt_dim2+ispace;
    v_log(cnt)=k_t(ispace);
        
    end
    
    end
end
x_log=x_log(1:cnt);
y_log=y_log(1:cnt);
v_log=v_log(1:cnt);

filt_mat=sparse(x_log,y_log,v_log,movie_time*1,filt_dim1*filt_dim2*movie_time,1*movie_time*filt_len*filt_dim1*filt_dim2);
inv_filt_mat=(filt_mat*filt_mat')\filt_mat;
end