function make_saccades_target_movie()
% MAKE_SACCADES_TARGET_MOVIE Generates a movie of saccade targets. The
% movie consists of a sequence of targets randomly displayed on the screen.
%   

% Set parameters
DISPLAY_FIGURE = false;
SAVE_MOVIE = true;
MOVIE_OUTPUT_PATH = '/Volumes/Lab/Projects/saccades-positioning/stimuli/jumping_target_int30.rawmovie';

% Set state of random number generator
par = saccades_movie_params();
s = RandStream('mt19937ar','Seed',par.RND_SEED);
RandStream.setGlobalStream(s);

% Generate target mask
target_mask = generate_saccade_target_mask();

% Generate sequence of target positions
tpos_seq = generate_sequence_target_pos(size(target_mask));

if SAVE_MOVIE
    % Generate raw movie header
    movie_length = (par.N_FRAMES_DISPLAY + par.N_FRAMES_INBETWEEN) * par.N_POS_TARGETS;
    rawmovie_header = sprintf('width\t%d\r\n', par.MOVIE_W);
    rawmovie_header = [rawmovie_header, sprintf('height\t%d\r\n', par.MOVIE_W)];
    rawmovie_header = [rawmovie_header, sprintf('frames-generated\t%d\r\n', movie_length)];
    rawmovie_header = [rawmovie_header, sprintf('interval\t%d\r\n', par.INTERVAL)];
    rawmovie_header = [rawmovie_header, sprintf('algorithm\tmake_saccades_target_movie\r\n')];
    rawmovie_header = [rawmovie_header, sprintf('ref\tvanBeers2007\r\n')];
    rawmovie_header = [rawmovie_header, sprintf('\r\n')];

    % Create rawmovie and write header to file
    fid = fopen(MOVIE_OUTPUT_PATH, 'w');  
    assert(fid ~= -1);
    fprintf(fid, sprintf('header-size\t%.10d\r\n', length(rawmovie_header)+24));
    fprintf(fid, rawmovie_header);
end
    
% Generate sequence of images
img_bgnd = ones(par.MOVIE_H, par.MOVIE_W)*par.BGND_COLOR;
if DISPLAY_FIGURE
    fh = figure();
    imshow(img_bgnd)
end
for kk=1:par.N_POS_TARGETS
    im_frame = img_bgnd;
    im_frame(tpos_seq(kk, 1)+(1:size(target_mask, 1)), ...
        tpos_seq(kk, 2)+(1:size(target_mask, 2))) = target_mask;
    
    for ll=1:par.N_FRAMES_DISPLAY
        if DISPLAY_FIGURE
            imshow(im_frame)
        end
        if SAVE_MOVIE
            % Convert image to 0-255 before saving
            im = round(reshape(im_frame,[1, numel(im_frame)])*par.MAX_GREY_LEVELS);
            fwrite(fid, repmat(im, [3,1]), 'ubit8');
        end
    end

    for ll=1:par.N_FRAMES_INBETWEEN
        if DISPLAY_FIGURE
            imshow(img_bgnd)
        end
        if SAVE_MOVIE
            % Convert image to 0-255 before saving
            im = round(reshape(img_bgnd,[1, numel(img_bgnd)])*par.MAX_GREY_LEVELS);
            fwrite(fid, repmat(im, [3,1]), 'ubit8');
        end
    end

end
    
if SAVE_MOVIE
    fclose(fid);
end

end % make_saccades_target_movie

function par = saccades_movie_params()
% SACCADES_MOVIE_PARAMS Sets the parameters for generating the saccade
% target movie.

par.MOVIE_W = 320;
par.MOVIE_H = 320;
par.MOVIE_STX_WIDTH = 1;
par.N_POS_TARGETS = 2500;
par.UM_TO_DEGREE = 200;
par.PIXEL_TO_UM = 4.5;
par.MAX_GREY_LEVELS = 255;

par.SZ_OUTER_RING_DEG = 0.33;
par.SZ_INNER_RING_DEG = 0.17;

par.BGND_COLOR = 0.98;
par.RING_COLOR = 0.02;

par.INTERVAL = 30;
par.N_FRAMES_DISPLAY = 3;
par.N_FRAMES_INBETWEEN = 1;

par.RND_SEED = 0;
par.UNIFORMITY_FACTOR = 1/5;

end % saccades_movie_params

function target_mask = generate_saccade_target_mask()
% GENERATE_SACCADE_TARGET_MASK Generates the shape of the saccade target.

par = saccades_movie_params();
mask_sz = ceil(par.SZ_OUTER_RING_DEG*par.UM_TO_DEGREE/par.PIXEL_TO_UM/2)*2 + 1;
outer_ring_rad = ceil(par.SZ_OUTER_RING_DEG*par.UM_TO_DEGREE/par.PIXEL_TO_UM/2);
inner_ring_rad = ceil(par.SZ_INNER_RING_DEG*par.UM_TO_DEGREE/par.PIXEL_TO_UM/2);

target_mask = par.BGND_COLOR*ones(mask_sz);

[xx, yy] = meshgrid((-(mask_sz-1)/2):1:((mask_sz-1)/2));
rr = sqrt(xx.^2 + yy.^2);
target_mask(rr < outer_ring_rad) = par.RING_COLOR;
target_mask(rr < inner_ring_rad) = par.BGND_COLOR;

end % generate_saccade_target_mask

function tpos_seq = generate_sequence_target_pos(tm_sz)
% GENERATE_SEQUENCE_TARGET_POS Generates a sequence of target positions for
% displaying the saccade target object.

par = saccades_movie_params();
max_y_pos = par.MOVIE_H - tm_sz(2);
max_x_pos = par.MOVIE_W - tm_sz(1);

% Position of the targets is biased towards the center of the array. 
% Make par.UNIFORMITY_FACTOR larger to get more uniform distributions
ypos_targets = ceil(randn(par.N_POS_TARGETS, 1)*max_y_pos*par.UNIFORMITY_FACTOR + max_y_pos/2);
ypos_targets= mod(ypos_targets, max_y_pos) + 1;

xpos_targets = ceil(randn(par.N_POS_TARGETS, 1)*max_x_pos*par.UNIFORMITY_FACTOR + max_x_pos/2);
xpos_targets = mod(xpos_targets, max_x_pos) + 1;

tpos_seq = [ypos_targets xpos_targets];
tpos_seq(tpos_seq == 0) = 1;

end % generate_sequence_target_pos