%% USED TO MAKE FIG 1
clear 
cell_type = 'On Parasol';
Corr_NS = [];
Loc = [];
Color_idx = [];
default_colors = get(gca,'ColorOrder');
stim_end = 20 - [240]/8+30;
%datarun{1} = load_data('/Volumes/Analysis/2015-09-23-0/data012-data013-data017-norefit/data017-from-data012_data013_data017/data017-from-data012_data013_data017');
%datarun{2} = load_data('/Volumes/Analysis/2015-09-23-0/data012-data013-data017-norefit/data012-from-data012_data013_data017/data012-from-data012_data013_data017');
%datarun{3} = load_data('/Volumes/Analysis/2015-09-23-0/data012-data013-data017-norefit/data013-from-data012_data013_data017/data013-from-data012_data013_data017');
stim_end = 20 - [280]/8+30;
datarun{1} = load_data('/Volumes/Analysis/2015-10-06-0/data000-data015-norefit/data000-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015/data000-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015');
datarun{2} = load_data('/Volumes/Analysis/2015-10-06-0/data000-data015-norefit/data004-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015/data004-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015');
datarun{3} = load_data('/Volumes/Analysis/2015-10-06-0/data000-data015-norefit/data005-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015/data005-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015');
% stim_end = 20 - [240,320]/8+30;
% datarun{1} = load_data('/Volumes/Analysis/2015-10-29-7/data008-data012/data012/data012');
% datarun{2} = load_data('/Volumes/Analysis/2015-10-29-7/data008-data012/data008/data008');
% datarun{3} = load_data('/Volumes/Analysis/2015-10-29-7/data008-data012/data009/data009');
% datarun{4} = load_data('/Volumes/Analysis/2015-10-29-7/data008-data012/data010/data010');
% datarun{1} = load_data('/Volumes/Analysis/2015-11-09-1/data009-data013/data009/data009');
% datarun{2} = load_data('/Volumes/Analysis/2015-11-09-1/data009-data013/data012/data012');
% datarun{3} = load_data('/Volumes/Analysis/2015-11-09-1/data009-data013/data013/data013');
%datarun{1} = load_data('/Volumes/Analysis/2015-10-29-2/data045-data048/data048/data048');
%datarun{2} = load_data('/Volumes/Analysis/2015-10-29-2/data045-data048/data045/data045');
%datarun{3} = load_data('/Volumes/Analysis/2015-10-29-2/data045-data048/data046/data046');
%datarun{4} = load_data('/Volumes/Analysis/2015-10-29-2/data045-data048/data047/data047');

datarun{1} = load_params(datarun{1});
datarun{1} = load_sta(datarun{1});
cids = get_cell_ids(datarun{1}, cell_type);
% cids = [cids get_cell_ids(datarun{1}, 'Off Parasol New Class')];
[NSEM_Corr,locations, avg_profile, PSTHs] = wide_field(datarun, cell_type,stim_end, cids);
[NSEM_Corr_AC,locations_AC,avg_profile, PSTHs_AC] = wide_field(datarun, 'Off Amacrine',stim_end, get_cell_ids(datarun{1}, 'Off Amacrine'));
datarun{1}.default_sta_fits = 'vision';
% [r,m,s] = average_radius(datarun{1}, cell_type);
Corr_NS = [Corr_NS; NSEM_Corr(:,3)];
%Loc = [Loc ; (locations(:,2)-stim_end)/(2*r)];
% AC_Loc = (locations_AC(:,2)-stim_end)/(2*r);
Color_idx = [Color_idx; ones(length(locations), 1)];

%%
figure(10); hold on; plot(locations(57:end), NSEM_Corr(:,4), '.', 'Color', default_colors(1,:), 'MarkerSize', 10)
title('MSE')

%%
default_colors = get(gca,'ColorOrder');
hold on
h = figure;
hold on
set(h, 'Position', [100 100 1000 2500]);
[~, idx] = sort(locations_AC);
for i = 1:3:length(NSEM_Corr_AC)
   plot(PSTHs_AC{idx(i),1}{3}/150+locations_AC(idx(i)), 'Color',default_colors(4,:))
end
% for i = 1:length(NSEM_Corr_AC)
%    plot(PSTHs_AC{i,1}{4}/1000+locations_AC(i+5), 'g')
% end
[~, idx] = sort(locations);
for i = 1:3:length(cids)
   %plot(PSTHs{idx(i),1}{2}/1000+locations(idx(i)),'Color', default_colors(1,:))
   plot(PSTHs{idx(i),1}{3}/150+locations(idx(i)),'Color', default_colors(2,:))
end
% for i = 1:length(cids)
%    plot(PSTHs{idx(i),1}{2}/1000+locations(idx(i)+56),'Color', default_colors(1,:))
%    plot(PSTHs{idx(i),1}{4}/1000+locations(idx(i)+56),'Color', default_colors(2,:))
% end
for i = 1:15
    plot([i, i]*100, [-8 8],'Color', default_colors(3,:))
end
% hold on
% plot(1550+avg_profile(:,1), avg_profile(:,2)-13, 'k')

ylim([-7 7])
xlim([-300 1800])
set(gca, 'XTick', [0 1500])
set(gca, 'XTickLabel', [0 15])

%%

cidx = [79, 95, 86];
h = figure;
hold on
set(h, 'Position', [100 100 600 200]);
for i = 1:15
   plot(i*100*[1 1],[0 200], 'Color', 0.8*[1 1 1], 'LineWidth', 2) 
end
for i = [cidx(1)]% 29 4]
   plot(PSTHs{i,1}{2},'k', 'LineWidth', 2)
   plot(PSTHs{i,1}{3},'Color', default_colors(2,:), 'LineWidth', 2)
end
xlim([0 500])
set(gca, 'XTick', [0 500]);
set(gca, 'XTickLabel', [0 5])
set(gca, 'YTick', [0 100])
ylim([0 150])

h = figure;
hold on
set(h, 'Position', [100 100 600 200]);
for i = 1:15
   plot(i*100*[1 1],[0 200], 'Color', 0.8*[1 1 1], 'LineWidth', 2) 
end
for i = cidx(2)
   plot(PSTHs{i,1}{2},'k', 'LineWidth', 2)
   plot(PSTHs{i,1}{3},'Color', default_colors(1,:), 'LineWidth', 2)
end

xlim([500 1000])
set(gca, 'XTick', [500 1000]);
set(gca, 'XTickLabel', [0 5])
set(gca, 'YTick', [0 100])
ylim([0 150])

h = figure;
hold on
set(h, 'Position', [100 100 600 200]);
for i = 1:15
   plot(i*100*[1 1],[0 200], 'Color', 0.8*[1 1 1], 'LineWidth', 2) 
end
for i = cidx(3)
   plot(PSTHs{i,1}{2},'k', 'LineWidth', 2)
   plot(PSTHs{i,1}{3},'Color', default_colors(5,:), 'LineWidth', 2)
end
xlim([1000 1500])
set(gca, 'XTick', [1000 1500]);
set(gca, 'XTickLabel', [0 5])
set(gca, 'YTick', [0 100])
ylim([0 150])

cids_examples = cids(cidx);

% Stimulus images
figure; 
load('/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian/matfiles/movie_chunk_1.mat');
temp = movie(:,:,1)';
imagesc(flipud(temp))
hold on
datarun = load_data('/Volumes/Analysis/2015-10-29-2/data045-data048/data048/data048');
datarun = load_params(datarun);
cids = get_cell_ids(datarun, cell_type);
cids = [cids get_cell_ids(datarun, 'Off Parasol New Class')];
plot_rf_fit(datarun, cids, 'scale', 4, 'fill', true, 'fill_color', [1 1 1], 'edge', false);
plot_rf_fit(datarun, [cids_examples(2)], 'scale', 4, 'fill', true, 'fill_color', default_colors(1,:), 'edge', false);
plot_rf_fit(datarun, [cids_examples(1)], 'scale', 4, 'fill', true, 'fill_color', default_colors(2,:), 'edge', false);
plot_rf_fit(datarun, [cids_examples(3)], 'scale', 4, 'fill', true, 'fill_color', default_colors(5,:), 'edge', false);
colormap gray
axis image
axis off

%hold on
%plot(325+avg_rf/5, 4*(1:200)/5, 'k', 'LineWidth', 2)

figure; 
load('/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian/matfiles/movie_chunk_1.mat');
temp(81:end,:) = 64;
imagesc(flipud(temp))
hold on
cids = get_cell_ids(datarun, cell_type);
cids = [cids get_cell_ids(datarun, 'Off Parasol New Class')];
plot_rf_fit(datarun, cids, 'scale', 4, 'fill', true, 'fill_color', [1 1 1], 'edge', false);
plot_rf_fit(datarun, [cids_examples(2)], 'scale', 4, 'fill', true, 'fill_color', default_colors(1,:), 'edge', false);
plot_rf_fit(datarun, [cids_examples(1)], 'scale', 4, 'fill', true, 'fill_color', default_colors(2,:), 'edge', false);
plot_rf_fit(datarun, [cids_examples(3)], 'scale', 4, 'fill', true, 'fill_color', default_colors(5,:), 'edge', false);
colormap gray
axis image
axis off



%{


%%
hold on
h = figure(1);
set(h, 'Position', [100 100 500 1000]);
for i = 1:109
   std_PSTH(i) = std(PSTHs{i,1}{3});
end
plot(locals, std_PSTH,'.')

%%

trials = 40;
seconds = 15;
bindur = 1/1000;

frames = seconds*120;
datarun{3} = load_neurons(datarun{3});
prepped_data = interleaved_data_prep(datarun{3}, frames, trials, 'cell_spec', cell_type, 'datarun_class', datarun{1});
bins =  ceil(seconds/bindur);
n_cells = length(cids);


%%
BPS = zeros(n_cells,1);
for i = 1:n_cells
    
    logicalspike = zeros(trials,bins) ;
    for i_blk = 1 : trials
        spt = prepped_data.testspikes{i_blk,i};
        binnumber = ceil(spt / bindur );
        logicalspike( i_blk, binnumber )  =  logicalspike( i_blk,binnumber ) + 1;
    end
    clear i_blk spt sptimes
    
    spikerate_bin    = size(find(logicalspike(:))) /  size(logicalspike(:));
    
    model_null0      = spikerate_bin * ones(1, bins);
    model_uop0       = (1/trials) * sum(logicalspike,1);
    model_null       = repmat(model_null0, trials, 1);
    model_uop        = repmat(model_uop0, trials, 1);
    
    null_logprob = sum(eval_rasterlogprob(logicalspike, model_null, 'binary', 'conditioned'));
    uop_logprob = sum(eval_rasterlogprob(logicalspike, model_uop, 'binary', 'conditioned'));
    
    uop_bits             = uop_logprob - null_logprob;
    BPS(i)    = uop_bits; % / (sum(model_null0));
    
end

%%
figure;
plot(locals, BPS, '.', 'MarkerSize', 10)
title('Bits')
%}