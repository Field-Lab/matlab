function val =dual_AAtx2(nu_vec,version)
global_vars

if(version==2) % A'x
nu=reshape(nu_vec,[movie_time,n_cell]);

val=Atx(stas,nu,filt_dim1,filt_dim2,filt_len,movie_time,n_cell);

val=val(:);
else % Ax
    nu=reshape(nu_vec,[filt_dim1,filt_dim2,movie_time]);
    val=Ax(stas,nu,movie_time,n_cell);
    val=val(:);

    
end
end