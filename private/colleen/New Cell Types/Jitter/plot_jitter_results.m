  number = [184 241 437 481 545 811 1386 1531 1786 1936 2011 2041 2103 2161 2331 2431 3048 3110 3166 3288 3436 3816 3886 3901 4113 4280 4550 4892 4937 4981 5012 5062 5071 5523 5915 6211 6331 6511 6526 7051 7456] %ON parasol
    number = [2011 2041 2103 2161 2331 2431 3048 3110 3166 3288 3436 3816 3886 3901 4113] %ON parasol
   number = [481 811 1531 2103 2161 3288 3816 3886 4280 4937 6331 6511] %ON parasol
   number = [185 1217 1221 1938 1942 3037 3904 3908 4893 4982 5523 5526 6094 6533 7056 7475 481 5061 3618 3811 5916 814 1477 1986 2167 2168 2991 3289 4326 4941 6513 6517 1441 1786 2431 3436 5026 6212 545 394 1218 1744 1970 2420 3302 4053 4097 4473 5987 6063 6617] %ON parasol
  number = [185 1217 1221 1938 1942 3037 3904 3908 4893 4982 5523 5526 6094 6533 7056 7475 481 5061 3618 3811 5916 814 1477 1986 2167 2168 2991 3289 4326 4941 6513 6517 1441 1786 2431 3436 5026 6212 545 394 1218 1744 1970 2420 3302 4053 4097 4473 5987 6063 6617] %ON parasol
number = [1711 3077 4337 4761 5888 3078 3381 4248 4653 5886 7507 7508];
number = [1279 2403 2584 2587 3380 5074 5104 5211 6128 6562 7326 7432 3124 4609 4684 4866 7028 7055 7221 1955 2990 4207 5852 6106 6391 6438 7187 396 992 4208 5194 6771 6968 7595 752 3318 4163];
number = [1279 2403 2584 2587 3380 5074 5104 5211 6128 6562 7326 7432 3124 4609 4684 4866 7028 7055 7221 1955 2990 4207 5852 6106 6391 6438 7187 396 992 4208 5194 6771 6968 7595 752 3318 4163];
number = [1021 1022 1023 1027 1028 1741 1742 1743 1744 2751 2752 2753 3066 3067 3068 3466 3468 3469 3470 4296 4297 4298 4299 4300 4327 4328 4329 4330 4332 5026 5027 5028 5176 5177 5178 5179 5180 5536 5537 5538 5539] %ON parasol
    number = [154 408 903 1953 3137 3259 3319 3721 3813 3861 4490 4565 4773 4778 5912 1711 3077 4337 4761 5888 3078 3381 4248 4653 5886 7507 7508] %ON parasol
 %% data026_cf_split/edited/data026_cf_split
    number = [1021 1022 1023 1027 1028 1741 1742 1743 1744 2751 2752 2753 3066 3067 3068 3466 3468 3469 3470 4296 4297 4298 4299 4300 4327 4328 4329 4330 4332 5026 5027 5028 5176 5177 5178 5179 5180 5536 5537 5538 5539] %ON parasol
   %ON parasol
    number = [2751 2752 2753, [nan], [nan];...
                3066 3067 3068, [nan], [nan];...
                3466 3468 3469 3470, [nan];...
                4296 4297 4298 4299 4300;...
                4327 4328 4329 4330 4332] 
%OFF parasol
 number = [1021 1022 1023 1027 1028;
           1741 1742 1743 1744 nan;
           5026 5027 5028 nan nan;
           5176 5177 5178 5179 5180;
           5536 5537 5538 5539 nan] 

          %%      
date = '2016-02-17-6';
datarun = 'data026_cf_split/data026_cf_split_original';
% number = [154 408 903 1953 3137 3259 3319 3721 3813 3861 4490 4565 4773 4778 5912];

% number = [532 1022 1568 2376 3186 4041 4703 5239 6788 7101];
% number{2} = [185 306 413 996 1058 1759 1821 2061 3396 3593 3830 4656 5342 5733 6064 6213 7148];

% number = [532 1022 1568 2376 3186 4041 4703 5239 6788 7101];
%     number= []; %ON parasol

%number = [1936] %ON parasol

  % dataparam.cell_specification = [482 813 1537 2103 2167 3288 3694 3889 4326 4939 6336 6517];

number = number';
sta = zeros(320,640);
figure;
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

x_plots = size(number,1); % number of plots down
y_plots = size(number,2);

% y_plots = ceil(length(number)/x_plots);
ha = tight_subplot(x_plots, y_plots, [.01 .03],[.01 .01],[.01 .01]);

for i = 1:x_plots*y_plots
    if ~isnan(number(i))%i <= length(number(:))
%     hold on 
    load(['/Volumes/Lab/Users/crhoades/Jitter/', date, '/', datarun,'/Cell ', num2str(number(i)), '.mat'])
    temp = permute(temp, [2,1,3,4]);
%     [sig_stixels] = significant_stixels(temp, 'select', 'thresh', 'thresh', 4.25);
    [~,start_index] = max(sum(reshape(temp.^2,[],size(temp,4)),1));

    temp = norm_image(temp);
    axes(ha(i));

    image(temp(:,:,:,start_index));
    axis image
    title(['Cell ', num2str(number(i))])
    axis off
    
%     sta((sig_stixels) == 1) = i;
%    sta = sta+ sig_stixels;

    else
            axes(ha(i));

        axis off
    end

    
end
datarun_mod = strrep(datarun, '_', '\_');
suptitle({[date, ' ', datarun_mod]; 'OFF Parasol'} )


figure; imagesc(sta(10:end-10, 10:end-10));

axis image
axis off


[~,git_hash_string] = system('git rev-parse HEAD');
fprintf('Git Hash: %s \n', git_hash_string);


%% Mosiac from Vision
dataparam.date='2015-09-23-7';
dataparam.concatname='data031';
dataparam.file_name = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
select_cells = 1;
if select_cells == 1
dataparam.cell_specification = [532 1022 1568 2376 3186 4041 4703 5239 6788 7101];
    
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
summary = zeros(size(double(datarun.stas.stas{cell_ids(1)}),1), size(double(datarun.stas.stas{cell_ids(1)}),2));

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

