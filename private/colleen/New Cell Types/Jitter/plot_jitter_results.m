  number = [184 241 437 481 545 811 1386 1531 1786 1936 2011 2041 2103 2161 2331 2431 3048 3110 3166 3288 3436 3816 3886 3901 4113 4280 4550 4892 4937 4981 5012 5062 5071 5523 5915 6211 6331 6511 6526 7051 7456] %ON parasol
    number = [2011 2041 2103 2161 2331 2431 3048 3110 3166 3288 3436 3816 3886 3901 4113] %ON parasol
   number = [481 811 1531 2103 2161 3288 3816 3886 4280 4937 6331 6511] %ON parasol
   number = [185 1217 1221 1938 1942 3037 3904 3908 4893 4982 5523 5526 6094 6533 7056 7475 481 5061 3618 3811 5916 814 1477 1986 2167 2168 2991 3289 4326 4941 6513 6517 1441 1786 2431 3436 5026 6212 545 394 1218 1744 1970 2420 3302 4053 4097 4473 5987 6063 6617] %ON parasol
  number = [185 1217 1221 1938 1942 3037 3904 3908 4893 4982 5523 5526 6094 6533 7056 7475 481 5061 3618 3811 5916 814 1477 1986 2167 2168 2991 3289 4326 4941 6513 6517 1441 1786 2431 3436 5026 6212 545 394 1218 1744 1970 2420 3302 4053 4097 4473 5987 6063 6617] %ON parasol

  %number = [1936] %ON parasol

  % dataparam.cell_specification = [482 813 1537 2103 2167 3288 3694 3889 4326 4939 6336 6517];


sta = zeros(320,640);
figure;
for i = 1:length(number)
    hold on 
    load(['/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data024/Cell ', num2str(number(i)), '.mat'])
    temp = permute(temp, [2,1,3,4]);
%     [sig_stixels] = significant_stixels(temp, 'select', 'thresh', 'thresh', 4.25);
     plot_sta_(temp)
    title({'2016-02-17-6 data024' ;['Cell ', num2str(number(i))]})
    axis off
%     sta((sig_stixels) == 1) = i;
%    sta = sta+ sig_stixels;

end


figure; imagesc(sta(10:end-10, 10:end-10));

axis image


[~,git_hash_string] = system('git rev-parse HEAD');
fprintf('Git Hash: %s \n', git_hash_string);


%% Mosiac from Vision
dataparam.date='2016-02-17-6';
dataparam.concatname='data025';
dataparam.file_name = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
select_cells = 1;
if select_cells == 1
dataparam.cell_specification = [482 813 1537 2103 2167 3288 3694 3889 4326 4939 6336 6517];
    
end
dataparam.cell_type = {'all'};

% Right Movie
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name, '.sta'];

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

cell_ids = get_cell_indices(datarun, dataparam.cell_specification)
clear sta
summary = zeros(20,40);

for i = 1:length(cell_ids)
sta = double(datarun.stas.stas{cell_ids(i)});

  [sig_stixels, params, rf_strength_out] = significant_stixels(sta, 'select', 'thresh', 'thresh', 4);

    summary = summary+ sig_stixels;

end
    figure; imagesc(summary);
    axis image
    
    [~,git_hash_string] = system('git rev-parse HEAD');
fprintf('Git Hash: %s \n', git_hash_string);


%% OFF smooth

number =[1936 1937 2014 3048 3903 3908 4113 4981 5013 5522  6110 6529 7054 7461 7462];
number =[1936 2011 3166 3901 4113 4981 5012 6526 7051 7456];

sta = zeros(320,640);
figure;
for i = 1:length(number)
    ref_size = zeros(size(sta));
    hold on 
    load(['/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data026_spikes/Cell ', num2str(number(i)), '.mat'])
    temp = permute(temp, [2,1,3,4]);
    [sig_stixels, params, rf_strength_out] = significant_stixels(temp, 'select', 'thresh', 'thresh', 3.5);

    sig_stixels = full(sig_stixels);
    sig_stixels = sig_stixels(10:end-10, 10:end-10);
biggestBlob_stixs = ExtractNLargestBlobs(sig_stixels, 1);

    ref_size(10:end-10, 10:end-10) = biggestBlob_stixs;
%           plot_sta_(temp)
%      title({'2016-02-17-6 data026' ;['Cell ', num2str(number(i))]})
%     axis off
    sta((ref_size) == 1) = i;
%     sta = sta+ ref_size;

end


figure; imagesc(sta(10:end-10, 10:end-10));

axis image


[~,git_hash_string] = system('git rev-parse HEAD');
fprintf('Git Hash: %s \n', git_hash_string);


%% Mosiac from Vision
dataparam.date='2016-02-17-6';
dataparam.concatname='data025';
dataparam.file_name = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
select_cells = 1;
if select_cells == 1
%     dataparam.cell_specification = [1936 1937  1938 1940 2012 3048 3903 3906 3907 3908 4114 4116 4981 4983 5013 5523 6111 6528 7055 7462];
number =[1936 1937 2014 3048 3903 3908 4113 4981 5013 5522  6110 6529 7054 7461 7462];

end
dataparam.cell_type = {'all'};

% Right Movie
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name, '.sta'];

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

cell_ids = get_cell_indices(datarun, dataparam.cell_specification)
clear sta
summary = zeros(20,40);

for i = 1:length(cell_ids)
sta = double(datarun.stas.stas{cell_ids(i)});

  [sig_stixels, params, rf_strength_out] = significant_stixels(sta, 'select', 'thresh', 'thresh', 4);

    summary = summary+ sig_stixels;

end
    figure; imagesc(summary);
    axis image
    
    [~,git_hash_string] = system('git rev-parse HEAD');
fprintf('Git Hash: %s \n', git_hash_string);

