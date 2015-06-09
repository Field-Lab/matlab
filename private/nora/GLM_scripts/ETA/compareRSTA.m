[U,S,V] = svd(reshape(eta_CP{15}, 80*40, 90));
[U2,S2,V2] = svd(reshape(eta_noCP{15}, 80*40, 90));
%%
s = 2;

%%
figure;
plot(V(:,s))
hold on
plot(V2(:,s))
hold off
legend('RSTA with CP', 'RSTA with noCP')
xlabel('Time')
title('Temporal Part of Filter (from SVD)')

%%
figure;
subplot(2,1,1)
imagesc(reshape(U(:,s), 80, 40)')
% caxis([-0.2 1])
axis image
colormap gray
subplot(2,1,2)
imagesc(reshape(U2(:,s), 80, 40)')
% caxis([-0.2 1])
axis image
colormap gray

