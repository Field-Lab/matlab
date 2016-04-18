%% compute correlation coefficient
width = 100;
height = 100;

stixel_size = 16;



sta = zeros(width, height);
P = zeros(width*height, width*height);
% pool = parpool(4);
for pixel = 2135% 1:width*height
    if mod(pixel, 100) == 0
        fprintf('%d out of %d \n', pixel, width*height);
    end
    
    [xref,yref] = ind2sub(size(sta), pixel);
    for i = 1:width*height
        [x,y] = ind2sub(size(sta), i);
        dx = abs(xref-x);
        dy = abs(yref-y);
        P(pixel, i) = max(((stixel_size-dx)/stixel_size),0)* max(((stixel_size-dy)/stixel_size),0);
    end
end


load('/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data026cf_spikes/Cell 1936.mat')
% Psparse = sparse(P);
% A = sqrtm(inv(Psparse));
% A = P^(-1/2);
temp_sta = temp(440:539, 200:299, 2,16);
white_sta = reshape(sqrtm(P)\temp_sta(:), 100, 100);
figure; imagesc(white_sta);

