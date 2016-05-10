% PIECE 1
% class_file = '2016-02-17-1/data019/data019';
% stixel_class = 8;
% mask_save_folder = '/Volumes/Lab/Users/Nora/new_stim_nora/mask_NSEM/2016-02-17/Maskin2/';
% mkdir(mask_save_folder);
% Make 2 groups with 3 on and 3 off parasols
% cells{1} = {7172, 4804, 2401};
% cells{2} = {5971, 1051, 3316};
% cells{1} = {7366, 4351, 2296};
% cells{2} = {6046, 1366, 3393};
% sigmas = [2 4 5 6];
% groupname = 'Maskin';
% groupname = 'Maskin2';

% % PIECE 2
% class_file = '2016-02-17-6/data000/data000';
% stixel_class = 8;
% mask_save_folder = '/Volumes/Lab/Users/Nora/new_stim_nora/mask_NSEM/2016-02-17-6/Maskin2/';
% mkdir(mask_save_folder);
% % Make 2 groups with 3 on and 3 off parasols
% cells{1} = {5929, 1113, 3668};
% cells{2} = {7414, 4592, 2342};
% % cells{1} = {6023, 1052, 3421};
% % cells{2} = {7522, 4653, 2221};
% sigmas = [2 4 5 6];
% % groupname = 'Maskin';
% groupname = 'Maskin2';

% % PIECE 8
% class_file = '2016-02-17-8/data000/data000';
% stixel_class = 8;
% groupname = 'Maskin_M';
% mask_save_folder = ['/Volumes/Lab/Users/Nora/new_stim_nora/mask_NSEM/2016-02-17-8/' groupname '/'];
% mkdir(mask_save_folder);

%Maskin
% cells{1} = {5973, 917, 3452};
% cells{2} = {7250, 4876, 2731};
% sigmas = 2 4 5 6

%MaskinA
% cells{1} = {4863};
% cells{2} = {2944};
% sigmas = [2 4 6 10];
% circular = 1, separate masks for each

% % PIECE 8
class_file = '2016-02-17-8/data006/data006';
stixel_class = 4;
groupname = 'Maskin_M';
mask_save_folder = ['/Volumes/Lab/Users/Nora/new_stim_nora/mask_NSEM/2016-02-17-8/' groupname '/'];
mkdir(mask_save_folder);
cells{1} = {1157, 1561, 2388, 3031, 3977, 4867, 5451, 5956, 6458, 7265};
sigmas = [2 4 6];

% calculate masks for all the sigmas. If the sigmas are too big, split into
% two groups
for i_sigma = sigmas
    disp(i_sigma)
    if i_sigma < 8
        mask = make_mask([cells{1}], class_file, 'sigmas', i_sigma, 'stixel_class', stixel_class);
        pause(0.001)
        save([mask_save_folder groupname '_allcells_sigma' num2str(i_sigma)], 'mask')
    else
        for i_cell = 1:length(cells)
            mask = make_mask(cells{i_cell}, class_file, 'sigmas', i_sigma, 'stixel_class', stixel_class);
            pause(0.001)
            save([mask_save_folder groupname '_cells' num2str(i_cell) '_sigma' num2str(i_sigma)], 'mask')
        end
    end
    
end
