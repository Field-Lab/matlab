% load image stack, compute max projection, save


switch 5
    case 1
        % image names
        im_path = '/snle/data/2007-09-18-0/confocal/2009-09-24/2007-09-18-0_63x_stitch_05_relevant/';
        name = @(z)sprintf('2007-09-18-0_63x_stitch_05_z%03d.tif',z);
        new_name = '2007-09-18-0_63x_stitch_05_max.tif';
        z_range = 2:30;
    case 2
        % image names
        im_path = '/snle/data/2007-09-18-6/confocal/2010-02-20/';
        name = @(z)sprintf('2007-09-18-6_63x_stitch_01_z%03d.tif',z);
        new_name = '2007-09-18-6_63x_stitch_01_max.tif';
        z_range = 5:45;
    case 3
        % image names
        im_path = '/snle/data/2007-09-18-6/confocal/2010-02-20b/';
        name = @(z)sprintf('2007-09-18-6_63x_stitch_02_z%03d.tif',z);
        new_name = '2007-09-18-6_63x_stitch_02_relevant_max.tif';
        z_range = 9:61;
    case 4
        % image names
        im_path = '/snle/data/2007-09-18-6/confocal/2010-02-21/';
        name = @(z)sprintf('2007-09-18-6_63x_stitch_03_z%03d.tif',z);
        new_name = '2007-09-18-6_63x_stitch_03_relevant_max.tif';
        z_range = 10:64;
    case 5
        % image names
        im_path = '/snle/data/2007-09-18-6/confocal/2010-02-21/';
        name = @(z)sprintf('2007-09-18-6_63x_stitch_04_z%03d.tif',z);
        new_name = '2007-09-18-6_63x_stitch_04_relevant_max.tif';
        z_range = 0:44;
end



% compute maximum


fprintf('Computing max for %d images from ''%s''\n%s\n',length(z_range),im_path,repmat('-',1,length(z_range)))
start_time_max_proj = clock;

% load first image
fprintf('.')
max_proj = imread([im_path name(z_range(1))]);
    
% load subsequent images
for zz=z_range(2:end)
    fprintf('.')
    % load source image
    max_proj = max(max_proj,imread([im_path name(zz)]));
end

fprintf(' done (%0.1f sec)\n\n',etime(clock,start_time_max_proj))

% save out max projection
fprintf('Saving final image as ''%s''...',new_name)
imwrite(max_proj,[im_path new_name],'tif')
fprintf('done\n\n')
