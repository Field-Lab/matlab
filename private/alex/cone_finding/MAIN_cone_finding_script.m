%% Input here

% main parameters
piece = '2012-09-13-2';
run = 'data001';
movie_descr = 'BW-2-5-0.48-11111-300x300-60.35.xml';
online = false;  % look in online folder (/Volumes/Acquisition/date)
streamed = true;  %  look for streamed data in offline folder (/Volumes/Analysis/date/streamed)

% additionbal parameters
nickname = '';
rig = 'A';
optical_path_direction = 'below';
display_type = 'oled';
cell_types = {1,2,3,4,5}; % to analyze for cone finding

%% load stuff

% load data
path2data = find_data(piece, run, streamed, online);
datarun = load_data(path2data);
% datarun = load_data('/Volumes/Analysis/2012-09-24-5/d03-06-07-norefit/data003/data003');
datarun = load_data('/Volumes/Analysis/2012-09-24-5/d00_06-norefit/data003-from-d00_06/data003-from-d00_06');
datarun = load_data('/Volumes/Analysis/2012-09-13-2/d01_09-norefit/data001-from-d01_09/data001-from-d01_09');

datarun.names.nickname = nickname;
datarun.piece.rig = rig;
datarun.piece.optical_path_direction = optical_path_direction;
datarun.piece.display = display_type;

datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);


% find movie, check if exists
movie_spec = fullfile(movies_path(), movie_descr);
if ~exist(movie_spec, 'file')
    fprintf('\nMOVIE DOESN''T EXIST\n')
end

% BW or RGB stimulus?
independent = strcmpi(datarun.stimulus.independent, 't');
field_width = datarun.stimulus.field_width;
field_height = datarun.stimulus.field_height;

info(datarun);

starttime = tic;

%% get sta info - summary and static nonlinearities

% robust_std_method is 1 to match old implementation.  Set to 3,5 for some speedup.
disp('Loading sta summaries and static NLs...')
datarun = get_sta_summaries(datarun, cell_types, ...
    'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct( ...
        'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
        'thresh',5,'robust_std_method',1));
    
% calculate static nonlinearities
datarun = load_java_movie(datarun, movie_spec);
start_time=0;
datarun = get_snls(datarun, datarun.cell_ids(get_cell_indices(datarun, cell_types)),'frames',-2:0,'start_time',start_time,'stimuli',10000,'new',true);



%% initial cone estimation, prior set up (local max)

% sig_stixel threshold in sigma
stixel_threshold = 4.5;

disp('Compute spatial sensitivity...');
combine_stixels = 'sum';
% combine_stixels = 'max';
[spatial_sensitivity,all_sig_stixels,spatial_cell_ids] = compute_spatial_sensitivity(datarun, cell_types, 'verbose', true, 'combine_stixels', combine_stixels, 'selection_params', struct('type', 'max', 'thresh', stixel_threshold));

thresh = 15;
thresh_spatial_sensitivity = spatial_sensitivity;
thresh_spatial_sensitivity(thresh_spatial_sensitivity < thresh) = 0;

surface = thresh_spatial_sensitivity;
radius = 1;

disp('Initial cone estimates...');

local_maxima = find_cones_in_rf(surface, 'filter', [], 'selection', struct('type', 'max', 'thresh', 0, 'radius', radius));
cones_labeled = bwlabel(local_maxima,8);
ncones = max(cones_labeled(:));

cone_centers = zeros(ncones,2);
for nn = 1:ncones
    [cone_centers(nn,1),cone_centers(nn,2)] = ait_centroid(cones_labeled == nn);
end
clear nn

figure
subplot(1,2,1)
imagesc(surface.^0.5)
colormap gray
hold on; plot(cone_centers(:,1), cone_centers(:,2), 'yo')

for i=1:length(cone_centers)
    
    D=pdist2(cone_centers(i,:),cone_centers);
    m(i)=min(D(D>0));
end

subplot(1,2,2)
hist(m,0:0.5:35)
axis([0 35 0 Inf])
title('Nearest Neighbor distance for cones. SET PRIOR!')

myPrior=[2.75,3.25]; % L-M cones prior

myPrior=[4.0,5.0]; % S cones prior

%% alternatively

datarun=load_sta(datarun,'load_sta','all');


kk=zeros(21,21,length(datarun.cell_ids));
cnt=1;
for i=1:length(datarun.cell_ids)
    tmp=sum(datarun.stas.stas{i}(:,:,:,6),3);
    [a,b]=max(abs(tmp(:)));
    [x,y]=ind2sub(size(tmp),b);
    if y>10 && x>10 && x<310 && y<310
        kk(:,:,cnt)=tmp(x-10:x+10,y-10:y+10)*sign(tmp(b));
        cnt=cnt+1;
    end
end
kk=kk(:,:,1:cnt-1);

figure
subplot(2,1,1)
tmp=mean(kk,3);
imagesc(tmp/max(tmp(:)))
subplot(2,1,2)
tmp=std(kk,0,3);
imagesc(tmp)

figure
tmp=mean(kk,3);
plot(tmp(1:end,11))
hold on
plot(tmp(11,1:end),'r')

myPrior=[4.0,4.25]; % L-M cones prior

%% set bcf params

bcf_params = struct;
        
switch independent
    case false
        % Fake a RGB STA from BW
        for cc = 1:length(datarun.stas.rfs)
            temp_rf = datarun.stas.rfs{cc};
            temp_rf = repmat(temp_rf, [1 1 3]);
            datarun.stas.rfs{cc} = temp_rf;
        end

        % BW use green
        bcf_params.kernel_plot_colors = ('g')';
        bcf_params.C_C = myPrior; % distance prior
        
        % kernel RGB
        bcf_params.kernel_colors = cone_rgb_expected(datarun);

    case true % RGB stimulus
        
        if cell_types{1}==5 % if only S cones:
                bcf_params.C_C = myPrior;
                bcf_params.kernel_plot_colors = ('b')';                
                for cc = get_cell_indices(datarun, cell_types)
                    temp_rf = datarun.stas.rfs{cc}(:,:,3);
                    temp_rf = repmat(temp_rf, [1 1 3]);
                    datarun.stas.rfs{cc} = temp_rf;
                end
                
                % kernel RGB
                bcf_params.kernel_colors.C = [1 1 1];
                independent=false;
                
        else            
        bcf_params.kernel_plot_colors = ('rgb')';
        bcf_params.LM_MM = myPrior;
        bcf_params.LM_S = [1.8 2.0];
        bcf_params.S_S = myPrior*1.5;
        
        % kernel RGB
        bcf_params.kernel_colors = cone_rgb_expected(datarun);
        
        end
end



% figure
bcf_params.plot_fig = 20;
bcf_params.cones_fig = 21;
bcf_params.dll_fig = 22;

% cell spec
bcf_params.cone_finding_cell_spec = cell_types;

% which cell types to identify the sampling of
bcf_params.regression_cell_spec = cell_types;

% size of subsets
bcf_params.padding_x = 5; % was 5 before
bcf_params.padding_y = 5; % was 5 before
% bcf_params.roi_x_size = 2*bcf_params.padding_x + 10; % was +10 before. in general, should be +(pad size)*2
% bcf_params.roi_y_size = 2*bcf_params.padding_y + 10; % was +10 before

bcf_params.roi_x_size = 20; % was +10 before. in general, should be +(pad size)*2
bcf_params.roi_y_size = 20; % was +10 before

% don't recompute if not necessary
bcf_params.new_W = 0;
bcf_params.new_STAs = 0;

% iterations per patch
bcf_params.num_iter = 100 ;

% radius of relevance
% sets which stixels (around the marks) are considered relevant in each STA
bcf_params.rel_radius = 4;



% start time of the datarun
bcf_params.start_time = start_time;

% arbitrary scale factor
bcf_params.magic_number = 1;

% cone density prior
bcf_params.q = 0.05;

% kernel spacing (in pixels) and radius
bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.65*[1 1 1];

% central square
% whole thing
bcf_params.relevant_region_x = [1 field_width-2*bcf_params.padding_x]; bcf_params.relevant_region_y = [1 field_height-2*bcf_params.padding_y];  % x == width, y == height??


% save('/Users/alexth/Desktop/variables_conefinding_bayes.mat','bcf_params','cell_indices','cell_types','datarun','extra_dirname_info',...
%    'field_height','field_width','independent','s_path')
%%
tic;
bcf = bayesian_cone_finding_loop(datarun,bcf_params);
toc

toc(starttime);


%%
% make some false cone locations and types if needed
datarun.cones.centers = [];
datarun.cones.types = [];
    
choose_magic_number(datarun,bcf,bcf_params);


%%
magic_number = 15
extra_dirname_info='';
save_bayesian_cones(datarun, bcf, bcf_params, magic_number, ['all_',extra_dirname_info], false, [], 'fit_foa', [], 'robust_std_method', 1);

%%

datarun = load_data('/Volumes/Analysis/2012-09-24-5/d00_06-norefit/data003-from-d00_06/data003-from-d00_06');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);
load('/Volumes/Analysis/2012-09-24-5/cone_data/d00_06-norefit_d-bayes-msf_15.00-all_/parameters.mat')
load('/Volumes/Analysis/2012-09-24-5/cone_data/d00_06-norefit_d-bayes-msf_15.00-all_/results.mat')

nickname = '';
rig = 'A';
optical_path_direction = 'below';
display_type = 'oled';
cell_types = {1,2,3,4,5}; % to analyze for cone finding

datarun.names.nickname = nickname;
datarun.piece.rig = rig;
datarun.piece.optical_path_direction = optical_path_direction;
datarun.piece.display = display_type;

%% Save out stuff for Jeremy Freeman's analysis
datarun = load_cones_ath(datarun,'d03'); % load_cones(datarun, 'Analysis');
datarun = load_neurons(datarun);
datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'd03');
conepreprocess_save(datarun, 'cone_data_ind', 'bayes');
