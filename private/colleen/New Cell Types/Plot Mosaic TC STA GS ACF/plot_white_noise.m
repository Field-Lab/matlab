% rand_sample = round(rand(16,16,3));

rand_sample = normrnd(0,0.5, [16,16,3])
upsamp_image(:,:,1)= imresize(rand_sample(:,:,1), 2, 'nearest');
upsamp_image(:,:,2) = imresize(rand_sample(:,:,2), 2, 'nearest');
upsamp_image(:,:,3) = imresize(rand_sample(:,:,3), 2, 'nearest');



figure; imagesc(upsamp_image)
axis equal
axis tight
axis off
