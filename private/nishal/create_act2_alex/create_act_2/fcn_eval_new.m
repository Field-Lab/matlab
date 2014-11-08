function val=fcn_eval_new(nu,stas,b)

filt_dim1=size(stas{1},1);
filt_dim2=size(stas{1},2);
filt_len=size(stas{1},4);
movie_time=size(nu,1);
n_cell=length(stas);

val=Atx(stas,nu,filt_dim1,filt_dim2,filt_len,movie_time,n_cell);
val=+0.25*(sum(sum(sum(val.^2))))+sum(sum(b.*nu));

end