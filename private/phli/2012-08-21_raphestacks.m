%% Setup
piece = '2012-08-21_OD_Raphe';

%%
% raphestacks{1} = build_stack('2012-11-16/2012-08-21 OD Raphe Tuj1/2012-08-21 OD Raphe Tuj1 stack.tif:[1-14]', 'basepath', server_data_path);
% save(fullfile(server_path(), piece, 'stacks'));

%%
load(fullfile(server_path(), piece, 'stacks'));

%% Simple Z-projections
raw = raw_stack(raphestacks{1});
raw = squeeze(raw);

maxproj = max(raw,[],3);
imwrite(maxproj, 'maxproj.tif');

meanproj = mean(raw,3);
normmeanproj = meanproj ./ max(meanproj(:));
imwrite(normmeanproj, 'meanproj.tif');

%% Reslicing
reslicestack = 1;
raphestacks{reslicestack} = load_slices(raphestacks{reslicestack});
stack_reslicer(raphestacks{reslicestack}, 'points_color', 'y', 'bgpoint_color', 'g', 'method', 'linear');
raphestacks{reslicestack}.reslice_points = gui_points;
save(fullfile(server_path(), piece, 'stacks'));