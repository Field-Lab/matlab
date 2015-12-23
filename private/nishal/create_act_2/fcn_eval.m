function val=fcn_eval(x,filt_mat,b)
val=(filt_mat'*x);
val=+0.25*(val'*val)+(b'*x);

end