function F = construct_fourier_basis(sidelen,maxf)


F = zeros(sidelen,sidelen,maxf^2);


A = fft(eye(sidelen));

k =1;
for i=1:maxf
    for j=1:maxf
        F(:,:,k) = A(:,i)*A(:,j)';
        k = k+1;
    end
end

F = reshape(F,sidelen^2,size(F,3));