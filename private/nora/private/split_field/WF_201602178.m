clear
Corr_NS = [];
Loc = [];
Color_idx = [];
modu=[];
no_modu= [];
default_colors = get(gca,'ColorOrder');
Analysis_Path = '/Volumes/Analysis/2016-02-17-8/';

% datarun 1 = class, 2 = full_rep, 3 = split
clear datarun NSEM_Corr locations stim_end
stim_end = [120, 240]/8+30;
dataruns{1} = load_data([Analysis_Path 'streamed/data000/data000']);
dataruns{2} = load_data([Analysis_Path 'map-from-data000/data001/data001']);

modu=[];
no_modu= [];

dataruns{1} = load_params(dataruns{1});
for i = 1:length(dataruns)
    dataruns{i} = load_neurons(dataruns{i});
end

example_cells = 0;
dataruns{1} = load_sta(dataruns{1});

%%

for cell_type = {'On Parasol'}%, 'Off Parasol'}
    
    % Plot mosaic
    if example_cells
        figure; plot_rf_fit(dataruns{1}, cell_type)
        ex_cells = get_cell_indices(dataruns{1}, example_cells);
        PSTH_ex_cells = cell(length(example_cells), length(dataruns)-1);
        for i=1:length(example_cells)
            hold on; plot_rf_fit(dataruns{1}, example_cells(i), 'fill', true)
        end
    end
    
    % avg RF
   
    avg_rf = sum(get_average_rf(dataruns{1}, cell_type, 'scale', 5),3);
    n_angles = 4;
    slice = zeros(201,1);
    range = -100:100;
    
    for i = 1:n_angles
        figure(1); imagesc(avg_rf); axis image; title(i)
        prof = improfile;
        prof = [zeros(100,1); prof; zeros(100,1)];
        [~,center] = max(prof);
        slice = slice + prof(center+range);
    end
    slice = slice/n_angles;
    plot(slice)
    
    a = fit((1:201)'/5,slice, fittype('gauss1'));
    width = a.c1/sqrt(2);
    avg_profile = [slice, (1:201)'/(5*width)];
    
    % Organize PSTH and calculate errors
    prepped_data = interleaved_data_prep(dataruns{2}, 1700, 120, 'cell_spec', cell_type{1}, 'datarun_class', dataruns{1}, 'visual_check', 0);
    cid = get_cell_indices(dataruns{1},cell_type{1});
    n_cells = size(prepped_data.testspikes,2);
    loc = zeros(n_cells,1);
    NSEM_Corr = zeros(n_cells, 2);
    i_ex_cell = 0;
    for i_cell = 1:n_cells
        idx = 1:40;
        for i = 1:3
            temp.testspikes = prepped_data.testspikes(idx,:);
            PSTH{i} = IDP_plot_PSTH(temp, i_cell);
            idx = idx+40;
        end
        for i_run = 2:3
            NSEM_Corr(i_cell,i_run) = err(PSTH{1}, PSTH{i_run})/1500;
        end
        loc(i_cell) = (dataruns{1}.vision.sta_fits{cid(i_cell)}.mean(2)-stim_end(1))/width;
        
        if any(example_cells) && any(ex_cells == cid(i_cell))
            i_ex_cell = i_ex_cell + 1;
            PSTH_ex_cells{i_ex_cell} = PSTH;
        end
    end
    
    %
    if length(stim_end) == 2
        loc = [loc; loc-diff(stim_end)/width];
    end
    
    Corr_NS = [Corr_NS; NSEM_Corr(:,2); NSEM_Corr(:,3)];
    Loc = [Loc ; loc];
    Color_idx = [Color_idx; 1*ones(length(loc),1)];
    
    
end

%%
cutoff = 50;
figure;  hold on;
%fill([prof1(:,2)-(95/8); prof1(:,2)-(95/8)]+4,0.01+[prof1(:,1)/10000; zeros(200,1)], 0.2+[0.7 0.7 0.7],'EdgeColor', 0.2+[0.7 0.7 0.7])
% x = avg_profile(:,2)-(95/8)+4;
% offset = 0.028;
% y = offset+avg_profile(:,1)/10000;
% plot(x(x>-1.5),y(x>-1.5), 'Color', [0.7 0.7 0.7]-0.2)
% plot(x(x>-1.5), zeros(length(x(x>-1.5)))+offset, 'Color', [0.7 0.7 0.7]-0.3);
plot(Loc(Loc>-cutoff), Corr_NS(Loc>-cutoff),'.', 'Color', 0.3*[1 1 1])
hold on;

plot([0 0], [-0.1 0.1], 'Color', 0.5*[1 1 1], 'LineWidth', 2)
[x, y, error, bincount, binedge] = curve_from_binning(Loc(Loc>-cutoff), Corr_NS(Loc>-cutoff), 'num_bins',8);
plot(x,y,'k', 'LineWidth', 2)

%ylim([-0.04 0.08])
set(gca, 'YTick', [0 0.03])
set(gca, 'XTick', [-5 0 5])
%xlim([-7.5 7.5])
%ylim([-0.001 0.039])

