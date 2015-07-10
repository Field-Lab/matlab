function [f_val,grad,Hess] = GMLM_quad_decomposed_fcn(x,mov_filtered,sta_f,N,dt,gamma)

k=x(1:end);

kkx = k'*mov_filtered;
ekx =(kkx'>0).*(kkx)'.^(gamma);

f_val = (sum(ekx)*dt/N) - k'*sta_f;

ekx2=gamma*(kkx'>0).*(kkx)'.^(gamma-1);

grad = mov_filtered*ekx2*dt/N - sta_f;

ekx3 = gamma*(gamma-1)*(kkx'>0).*(kkx)'.^(gamma-2);
Hess = (repmat(ekx3',[size(mov_filtered,1),1]).*mov_filtered)*mov_filtered'*(dt/N);

f_val=gather(f_val);
grad = gather(grad);
Hess = gather(Hess);
end
