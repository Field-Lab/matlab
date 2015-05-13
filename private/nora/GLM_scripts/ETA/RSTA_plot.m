s = 2;
subplot(2,1,1)
title('SVD of RSTA')
imagesc(reshape(U(:,s), 80, 40)'); axis image; colormap gray
subplot(2,1,2)
plot(V(:,s))