%% load data
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/Classification/
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/DS' cell analysis'/
cd ~/Desktop/DS_code/
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
% dg := drifting gratings
datadg{1} = load_data('/Volumes/Vivaldi/Analysis/2015-09-01-0/data005/data005', opt);
datadg{1}.names.stimulus_path = '/Volumes/Vivaldi/Analysis/2015-09-01-0/stimuli/s05.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);

datadg{2} = load_data('/Volumes/Vivaldi/Analysis/2015-09-01-0/data006/data006', opt);
datadg{2}.names.stimulus_path = '/Volumes/Vivaldi/Analysis/2015-09-01-0/stimuli/s06.mat';
datadg{2} = load_stim_matlab(datadg{2}, 'user_defined_trigger_interval', 10);

datadg{3} = load_data('/Volumes/Vivaldi/Analysis/2015-09-01-0/data007/data007', opt);
datadg{3}.names.stimulus_path = '/Volumes/Vivaldi/Analysis/2015-09-01-0/stimuli/s07.mat';
datadg{3} = load_stim_matlab(datadg{3}, 'user_defined_trigger_interval', 10);

datadg{4} = load_data('/Volumes/Vivaldi/Analysis/2015-09-01-0/data008/data008', opt);
datadg{4}.names.stimulus_path = '/Volumes/Vivaldi/Analysis/2015-09-01-0/stimuli/s08.mat';
datadg{4} = load_stim_matlab(datadg{4}, 'user_defined_trigger_interval', 10);

% mb := moving bars

datamb{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-09-01-0/data001/data001', opt);
datamb{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-09-01-0/stimuli/s01.mat';
datamb{1} = load_stim_matlab(datamb{1});

datamb{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-09-01-0/data002/data002', opt);
datamb{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-09-01-0/stimuli/s02.mat';
datamb{2} = load_stim_matlab(datamb{2});

datamb{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-09-01-0/data003/data003', opt);
datamb{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-09-01-0/stimuli/s03.mat';
datamb{3} = load_stim_matlab(datamb{3});

datamb{4} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-09-01-0/data004/data004', opt);
datamb{4}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-09-01-0/stimuli/s04.mat';
datamb{4} = load_stim_matlab(datamb{4});

% ffp := full field pulses
dataffp{1} = load_data('/Volumes/Vivaldi/Analysis/2015-06-09-0/data010/data010', opt);
dataffp{1}.triggers = dataffp{1}.triggers(2:end);


%%
temp_datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-06-09-0/data001/data001', opt);
temp_datarun.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-06-09-0/stimuli/s01.mat';
temp_datarun = load_stim_matlab(temp_datarun);

%% classify DS vs non-DS cells (using moving bars)
  
    i = 1; % which datarun to use for classification
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb_gdf(datamb{i},datamb{i}.cell_ids, datamb{1}.triggers(end)+6, 1);
    ds_struct = mbcellanalysis(NumSpikesCell, StimComb);

    a = 1; b = 2; % which parameters to use for classification
    Taylor1 = ds_struct.MAG{a,1}./((sum(ds_struct.RHO{a,1},2))');
    Taylor2 = ds_struct.MAG{b,1}./((sum(ds_struct.RHO{b,1},2))');
    
    figure
    plot(log(Taylor1), log(Taylor2), 'o')
    
    title('ndf0')
    xlabel('TP 1')
    ylabel('TP 2')
    hold on
    
    [x, y] = ginput;
    plot(x, y);
    IN = inpolygon(log(Taylor1), log(Taylor2), x, y);
    [~, I] = find(IN == 1);
    id_init = datamb{i}.cell_ids(I);

    [C ia ib] = intersect(id_init, datamb{i}.cell_ids);
    vc = ones(length(datamb{i}.cell_ids),1);
    vc(ib) = 2; %initializing ds cells to cluster 2, everything else cluster 1

    close all;
    X = [];
    N = [];
    p = [];
    X(:,1) = log(Taylor1)';
    X(:,2) = log(Taylor2)';
    [idx N p] = clustering_analysis_plots(X, 0,1, 2, 0, 1, 0, 0, 0,0, vc);

    ds_id = [];
    ds_id = datamb{i}.cell_ids(idx==2);
    nonds_id = datamb{i}.cell_ids(idx==1);
    
%% drifting grating

n = 1; % number of datarun

[raster_dg, DG, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, StimComb] = get_spikescellstim(datadg{i},ds_id,0);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster_dg{i} = get_ds_raster(datadg{i}, ds_id);
end

delta_p = 1; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
    MAG_all_norm_dg{i} = normalize_MAG(DG{i});
    rep = datadg{i}.stimulus.repetitions;
end

ll = {'NDF0'};

%% plot cell summary

for cc = 1:2 %length(ds_id)
    plot_ds_raster_RS(DG, raster_dg, cc, ds_id(cc), ll, 1, 1, 0)
end

%% classify DSGC into subtypes (directions)
d = 3; % which datarun to use for classification
t = 2; % which parameter to use for classification
h = figure;
dirn = 4; % number of preferred directions
set(h, 'Position', [1 1 1080 500])
compass(MB{d}.U{t}, MB{d}.V{t});
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(MB{d}.U{t}, MB{d}.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id(idx_dir{i});
end

%% DS tuning curves (drifting grating)
% all ds cells

dirn = 4; %number of preferred direction
D = 1; % which datarun to use to calculate the preferred direction

for d = 1:n
    p_direction = DG{D}.angle{delta_p}';
    xx = 0:pi/6:11*pi/6; % if number of direction is not 8, change the code here!!!
    xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 12); % and here!!!
    xx(xx>pi) = xx(xx>pi)-2*pi;
    xx(xx<-pi) = xx(xx<-pi)+2*pi;


    subplot(1, 1, d)
    for i = 1:dirn
        for cc = 1:length(ds_id)
            [xsort, seq] = sort(xx(cc, :));
            y_temp = DG{d}.rho{delta_p}(cc, :);
            plot(xsort, y_temp(seq), 'b')
            hold on
        end
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{d})
    xlim([-pi pi])
end

%subtypes
clear rho_dg_mean rho_dg_ste dsi_dg_mean dsi_dg_ste
for d = 1:n
    subplot(1, 1, d)
    for i = 1:dirn
        rho_dg{d}{i} = [];
        dsi_dg{d}{i} = [];
        for cc = 1:length(idx_dir{i})
            %if ~dg_idx(idx_dir{i}(cc), d) && sum(DG{d}.rho{delta_p}(idx_dir{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
            y_temp = DG{d}.rho{delta_p}(idx_dir{i}(cc), :);
            plot(xsort, y_temp(seq), color(i))
            ylim([0 1])
%             pause
            hold on
            rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
            dsi_dg{d}{i} = [dsi_dg{d}{i}; DG{d}.dsindex{delta_p}(idx_dir{i}(cc))];
            %end
        end
        rho_dg_mean{d}(i, :) = mean(rho_dg{d}{i});
        rho_dg_ste{d}(i, :) = std(rho_dg{d}{i})/sqrt(size(rho_dg{d}{i}, 1));
        dsi_dg_mean{d}(i) = mean(dsi_dg{d}{i});
        dsi_dg_ste{d}(i) = std(dsi_dg{d}{i})/sqrt(length(dsi_dg{d}{i}));
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{d})
    xlim([-pi pi])
end
dsi_dg_mean = cell2mat(dsi_dg_mean');
dsi_dg_ste = cell2mat(dsi_dg_ste');
% plot average (cell type)
figure
for d = 1:1
    subplot(1, 1, d)
    for i = 1:4
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(i));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ll{d});
end
ct = {'superior','nasal','inferior','temporal'};
legend(ct)

% plot average (light level)
figure
for i = 1:4
    subplot(2, 2, i)
    for d = 1:1
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(d));
        hold on
        % fwhm.m from Patrick Egan Matlab file exchange
        width(i) = fwhm(xsort,rho_dg_mean{d}(i,:));
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ct{i})
end
legend(ll)
% DSI
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = {'temporal'; 'inferior'; 'nasal'; 'superior'};
model_series = [dsi_dg_mean(1,1); dsi_dg_mean(1,2); dsi_dg_mean(1,3); dsi_dg_mean(1,4)];   
model_error = [dsi_dg_ste(1,1); dsi_dg_ste(1,2); dsi_dg_ste(1,3); dsi_dg_ste(1,4)];
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('DSI')
legend( 'NDF0');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end


%% Moving bars (modified to look at the 3rd data run of the MB stimulus only)
n=1;
i =1;
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx] = deal(cell(n, 1));
%for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb_gdf(datamb{i},ds_id,0,1);
    k = 1;
    MB{k} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb));
    raster_mb{k} = get_mb_raster(datamb{i}, ds_id, duration);
    %for j = 1:length(raster_mb{i})
        %if(mb_idx(j, i))
            %raster_mb{i}{j} = [];
        %end
    %end
    trial_dur{k} = get_mb_trial_dur(datamb{i});
%end

delta_p = 1; % choose which params to use to calculate preferred direction indices 
D = 1; % choose which datarun to use to calculate preferred direction indices
MAG_all_norm_mb = cell(n, 1);

for i = 1:n
    [raster_p_sum_mb{i}, p_idx{i}] = get_pdirection_raster(raster_mb{i}, MB{D}.angle{delta_p});
    MAG_all_norm_mb{i} = normalize_MAG(MB{i});
    rep = datamb{3}.stimulus.repetitions;
end

ll = {'NDF0'};
ct = {'temporal', 'inferior', 'nasal', 'superior'};

%% plot cell summary

for d = 1:1
    for cc = 6:7 %length(ds_id)
        plot_mb_raster_RS(MB, raster_mb, trial_dur, idx_dir{d}(cc), ds_id(idx_dir{d}(cc)), ll, 1, 1, 0)
    end
end

%% DS tuning curves (moving bar)
% all ds cells

dirn = 4;
D = 1;
T = 2;

for d = 1:1
    p_direction = MB{D}.angle{T}';
    xx = 0:pi/6:11*pi/6; % if number of direction is not 8, change the code here!!!
    xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 12); % and here!!!
    xx(xx>pi) = xx(xx>pi)-2*pi;
    xx(xx<-pi) = xx(xx<-pi)+2*pi;


   
    subplot(1, 1, d)
    for i = 1:dirn
        for cc = 1:length(ds_id)
            %if ~mb_idx(cc, d)
            [xsort, seq] = sort(xx(cc, :));
            y_temp = MB{d}.rho{T}(cc, :);
            plot(xsort, y_temp(seq), 'b')
            hold on
            %end
        end
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{d})
    xlim([-pi pi])
end

%subtypes
clear rho_mb_mean rho_mb_ste dsi_mb_mean dsi_mb_ste
for d = 1:1
    subplot(1, 1, d)
    for i = 1:dirn
        rho_mb{d}{i} = [];
        dsi_mb{d}{i} = [];
        for cc = 1:length(idx_dir{i})
            [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
            y_temp = MB{d}.rho{T}(idx_dir{i}(cc), :);
            plot(xsort, y_temp(seq), color(i))
            ylim([0 1])
            hold on
            rho_mb{d}{i} = [rho_mb{d}{i}; y_temp(seq)];
            dsi_mb{d}{i} = [dsi_mb{d}{i}; MB{d}.dsindex{T}(idx_dir{i}(cc))];
        end
        rho_mb_mean{d}(i, :) = mean(rho_mb{d}{i});
        rho_mb_ste{d}(i, :) = std(rho_mb{d}{i})/sqrt(size(rho_mb{d}{i}, 1));
        dsi_mb_mean{d}(i) = mean(dsi_mb{d}{i});
        dsi_mb_ste{d}(i) = std(dsi_mb{d}{i})/sqrt(length(dsi_mb{d}{i}));
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{d})
    xlim([-pi pi])
end
dsi_mb_mean = cell2mat(dsi_mb_mean');
dsi_mb_ste = cell2mat(dsi_mb_ste');
% plot average (cell type)
figure
for d = 1:1
    subplot(1, 1, d)
    for i = 1:4
        errorbar(xsort/pi*180, rho_mb_mean{d}(i, :), rho_mb_ste{d}(i, :), color(i));
        hold on
    end
    xlabel('degrees')
    ylabel('normalized mean spike number')
    title(ll{d});
end
legend(ct)

% plot average (light level)
figure
for i = 1:4
    subplot(2, 2, i)
    for d = 1:1
        errorbar(xsort/pi*180, rho_mb_mean{d}(i, :), rho_mb_ste{d}(i, :), color(d));
        hold on
    end
    ylim([0 1])
    xlabel('degrees')
    ylabel('normalized average response')
    title(ct{i})
end
legend(ll)
% DSI
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = {'unknown'; 'inferior'; 'unknown'; 'superior'};
model_series = [dsi_mb_mean(1,1); dsi_mb_mean(1,2); dsi_mb_mean(1,3); dsi_mb_mean(1,4)];   
model_error = [dsi_mb_ste(1,1); dsi_mb_ste(1,2); dsi_mb_ste(1,3); dsi_mb_ste(1,4)];
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('DSI')
legend( 'NDF0');
hold on;

numgroups = size(model_series, 1);  
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));


for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

%% full field pulses

n_ffp = 1; % number of datarun

[raster_ff, raster_ff_all] = deal(cell(n_ffp, 1));
for d = 1:n_ffp
    [raster_ff{d}, raster_ff_all{d}] = get_ffp_raster(dataffp{d}, ds_id, 3);
end

for i = 1:4 %length(ds_id) 
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 1800 800])
        for d = 1:n_ffp
            subplot(1, n_ffp, d)
            plot_ffp(raster_ff{d}, raster_ff_all{d}, i, 3)
            title([num2str(ds_id(i)) ' ' ll{d}])
        end
        
%         print_close(1, [24, 12], num2str(ds_id(i)))
end
