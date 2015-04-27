% Given an M X N matrix A we want to block A into an M x ceil(N/k) matrix by
% summing across blocks of k columns.

function B = block_spikes(A,k)

A = full(A);

[M N] = size(A);

Nbyk = N/k;


%fprintf('M=%d,N=%d,k=%d\n',M,N,k);

Bhat = reshape(A(:,1:floor(Nbyk)*k),M,k,floor(Nbyk));
Bhat = sum(Bhat,2);
B = reshape(Bhat,M,floor(Nbyk));
1;
if (~isempty(k*floor(Nbyk)+1:size(A,2)))
    B = [B sum(A(:,k*floor(Nbyk)+1:end),2)];
end



