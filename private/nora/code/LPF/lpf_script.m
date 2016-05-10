% where to save the matfiles
matfile_save_location = '/Volumes/Lab/Users/Nora/LPF/matfiles/';

disp('Making matfiles')
% make matfile chunks
[~, p] = lpf_make_mat(matfile_save_location, 'seconds', 3600, 'sigma', 2, 'tau', 6, 'back_rgb', 0.25, 'intensity_SD', 0.08);

disp('Making rawmovie')
% get params for the movie name
sd = num2str(p.Results.intensity_SD);
movie_path = ['/Volumes/Lab/Users/Nora/LPF/rawmovie/lpf-' num2str(p.Results.sigma) '-' num2str(p.Results.tau) '-0_' sd(3:4) '.rawMovie'];
clear sd

% make rawmovie
lpf_make_rawmovie(matfile_save_location, movie_path, p.Results.seconds, p.Results.dim);