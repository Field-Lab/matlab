function val =dual_AAtx(nu_vec,version)
global_vars


nu=reshape(nu_vec,[movie_time,n_cell]);

val=Ax(stas,Atx(stas,nu,filt_dim1,filt_dim2,filt_len,movie_time,n_cell),movie_time,n_cell);

val=val(:);

end