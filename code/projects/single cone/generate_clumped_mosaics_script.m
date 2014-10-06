% get the paths to cones.txt for all 12 datasets
% the get_dataset_folder script follows the fruit name links and makes sure
% that 'cones.txt' exists in all folders before returning a value

%[path, name] = get_dataset_folder('all','suffix','cones.txt');
%[path, name] = get_dataset_folder('plantain','suffix','cones.txt');

[LMS_paths, LMS_names, cone_paths] = get_LMS_paths('high', 'cone_finding', 'standard');

num_sets = length(cone_paths);
for ds = 1:num_sets
    temp_path = cone_paths{ds};
    temp_path = ['/marte/snle/lab/Experiments/Array/Shared/one/', temp_path, '/cones.txt'];
    path{ds} = temp_path;
end
name = LMS_names;

path_to_save = '/snle/lab/Experiments/Array/Shared/one/simulations';

% how many text files to generate of each type? (should be 64)
nMosaics = 16;

for n = 1:length(name)
    
    % permute 100 times and compare to observed clumping index
    % generate_clumped_mosaic_using_ci(path{n},name{n},'num_mosaics',1,'fit_permuted',false,'make_plots',true,'iterations',100);

    % save out permuted mosaics. do not fit: just permute once and save
%     generate_clumped_mosaic_using_ci(path{n},name{n},'save_permuted_only',true,...
%         'num_mosaics',nMosaics,'folder','permuted','save_text',true,...
%         'base_path', path_to_save);
    
    % fit clumping index of permuted mosaic to that of the observed data; generate 64 text files for each dataset
    generate_clumped_mosaic_using_ci(path{n},name{n},'fit_permuted',true,...
        'num_mosaics',nMosaics,'folder','extra-clumped','save_text',true,...
        'base_path', path_to_save, 'clumping_multiplier', 4.5);
    
    % fit clumping index of permuted mosaic to that of the observed data; swap neighbor cones only; generate 64 text files for each dataset
%    generate_clumped_mosaic_using_ci(path{n},name{n},'neighbor_clumped',true,...
%        'num_mosaics',nMosaics,'fit_permuted',true,'folder','neighbor','save_text',true);
    
    % fit drps of permuted mosaic to those of the observed data; generate 64 text files for each dataset
%    generate_clumped_mosaic_using_drp(path{n},name{n},...
%        'num_mosaics',nMosaics,'folder','drp','save_text',true,'make_plots',true);
    
end


