% generate stimulus

image = 0.5*ones(400,600,3);


image(:, 200:231,:) =  repmat(1, 400,32,3);



G = fspecial('Gaussian', [10 10],10);
Ig = imfilter(image,G,'same');


imageToSave = Ig(20:380, 20:580, 3);
% figure; imshow(imageToSave)
figure; imshow(imageToSave)
annotation('arrow',[0.37 0.45], [0.55 0.55], 'Color', 'k', 'LineWidth', 2)


set(gcf,'color','w');

axis image
save('stimulusDepiction', 'imageToSave')