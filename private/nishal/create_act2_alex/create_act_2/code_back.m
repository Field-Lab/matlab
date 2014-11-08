% making sparse toeplitz ?
%filt_mat= sparse([],[],[],movie_time*n_cell,filt_dim1*filt_dim2*movie_time,n_cell*movie_time*6*320*160);
%filt_mat= sparse([],[],[],movie_time*1,filt_dim1*filt_dim2*movie_time,n_cell*movie_time*6*320*160);
filt_mat=cell(n_cell,1);

for icell=1:n_cell
k=stas{icell};
nzero=1*movie_time*6*320*160;
x_log=zeros(nzero,1);
y_log=zeros(nzero,1);
v_log=zeros(nzero,1);
cnt=0;

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
x_log=x_log(1:cnt);
y_log=y_log(1:cnt);
v_log=v_log(1:cnt);

filt_mat{icell}=sparse(x_log,y_log,v_log,movie_time*n_cell,filt_dim1*filt_dim2*movie_time,n_cell*movie_time*6*320*160);
end

%%

mov_modify=mov_reshape;
% Alternating projections same as Dykstra's for multiple affine sub-spaces!
%mov_modify_sum=zeros(length(mov_reshape),1);
residual=zeros(n_cell,1);
for outer_iter=1:1
for icell=1:n_cell
    [outer_iter,icell]
    mov_modify=mov_modify- filt_mat{icell}'*(inv_filt_mat{icell}*mov_modify);
end

for icell=1:n_cell
residual(icell)=norm(filt_mat{icell}*mov_modify);
end

end

orig=zeros(n_cell,1);
for icell=1:n_cell
    orig(icell)=norm(filt_mat{icell}*mov_reshape);
end
% If constraints inconsistent, then least square? pinv? pseudo-norm ? 
