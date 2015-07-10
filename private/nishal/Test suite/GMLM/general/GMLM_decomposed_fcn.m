function [f_val,grad,Hess] = GMLM_decomposed_fcn(x,mov_filtered,sta_f,N,dt)

k=x(1:end);


ekx = exp(k'*mov_filtered )';
f_val = (sum(ekx)*dt/N) - k'*sta_f;

grad = mov_filtered*ekx *dt/N - sta_f;

Hess = (repmat(ekx',[size(mov_filtered,1),1]).*mov_filtered)*mov_filtered'*(dt/N);

f_val=gather(f_val);
grad = gather(grad);
Hess = gather(Hess);
end
