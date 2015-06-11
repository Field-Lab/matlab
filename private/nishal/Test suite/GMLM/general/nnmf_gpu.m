function [W,H] = nnmf_gpu(V,k,tol)

m=size(V,1);
n=size(V,2);

V=gpuArray(V);
W=gpuArray(rand(m,k)); H =gpuArray(rand(k,n));

error = norm(V-W*H,'fro')/norm(V,'fro');
while error>tol
H = H .* (W'*V)./((W'*W)*H);
W = W .* (V*H')./(W*(H*H'));

%H=Hnew; W=Wnew;
error = norm(V-W*H,'fro')/norm(V,'fro')
end

W=gather(W);
H=gather(H);

end