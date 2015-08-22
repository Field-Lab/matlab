function my_movie = compute_raw_wn_movie(rgb_flag, bin_flag, rgb, seed, width, height, probability, frames)

rng_init = Init_RNG_JavaStyle(seed);

% noise type
% 0 - binary BW; n_bits = 1
% 1 - binary RGB; n_bits = 3
% 2 - gaussian BW; n_bits = 8 (1 draw)
% 3 - gaussian RGB; n_bits = 8 (3 draws)
if strcmp(bin_flag, 'bin') % binary
    if strcmp(rgb_flag, 'rgb') % RGB
        noise_type = 1;
        n_bits = 3;
        tmp = [ 1 1 1;  1 1 -1;  1 -1 1;  1 -1 -1;...
            -1 1 1;  -1 1 -1;  -1 -1 1;  -1 -1 -1];
        tmp = tmp .* rgb + 0.5;
    else % BW
        noise_type = 0;
        n_bits = 1;
        tmp = [1 1 1; -1 -1 -1] * rgb + 0.5;
    end
else % gaussian
    if strcmp(rgb_flag, 'rgb')  % RGB
        noise_type = 3;
    else % BW
        noise_type = 2;
    end
    n_bits = 8;
    tmp = norminv((1:256)/257, 0, 1)';
    tmp = repmat(tmp,1,3) .* rgb + 0.5;
end

tmp = uint8(round(255 * tmp))';
lut = tmp(:);
map_back_rgb = uint8( round( 255 * [0.5 0.5 0.5]'));
map =[];

my_movie = zeros(width, height, 3, frames);

for i=1:frames
    img_frame = Draw_Random_Frame_opt(rng_init, width, height, lut,...
        map, map_back_rgb, width, height, noise_type, n_bits, probability);
    my_movie(:,:,:,i) = shiftdim(img_frame(1:3,:,:),1);
    
end

