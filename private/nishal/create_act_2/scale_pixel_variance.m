function [x_k_1,ratio] = scale_pixel_variance(x_k_1,mov_show)

cMap_mov_show  = sqrt(sum(mov_show.^2,3))/size(mov_show,3);
cMap_xk1 = sqrt(sum(x_k_1.^2,3))/size(x_k_1,3);

ratio = cMap_mov_show./cMap_xk1;
x_k_1 = x_k_1 .*repmat(ratio,[1,1,size(x_k_1,3)]);
cMap_xk1_corrected = sqrt(sum(x_k_1.^2,3))/size(x_k_1,3);

figure;
subplot(2,2,1);
imagesc(cMap_mov_show);
colorbar

subplot(2,2,2);
imagesc(cMap_xk1);
colorbar

subplot(2,2,3);
imagesc(cMap_xk1_corrected);
colorbar
end