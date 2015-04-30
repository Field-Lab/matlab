clear

datapath = '/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/fixedSP_rk1_linear_MU_PS_CP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{2}='NSEM_mapPRJ/';
fittypepath{1}='WN_mapPRJ/';

cell_count = 0;
total_cell_count = 0;

% coupling.info = cellid, exp, BPS WN, BPS NSEM
% coupling.CP = paired cell, CP_WN, CP_NSEM

for exp = 2
    for fittype = 1:2
        cell_count = total_cell_count;
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) '*.mat']);
        n_cells=length(matfiles);
        
        disp(['Experiment ' num2str(exp) ' Stimulus ' fittypepath{fittype} ' with ' num2str(n_cells) ' Cells'])

        % Load cells
        for file=1:n_cells
            cell_count = cell_count + 1;
            load([datapath fittypepath{fittype} exp_names(exp,:) matfiles(file).name]);
            coupling{cell_count}.info(1) = str2double(matfiles(file).name(7:(end-4)));
            coupling{cell_count}.info(2+fittype) = fittedGLM.xvalperformance.glm_normedbits;
            coupling{cell_count}.CP(:,1) = double(fittedGLM.cellinfo.pairs);
            
            % Go through the pairs and extract coupling
            for pair = 1:6
                coupling{cell_count}.CP(pair,1+fittype) = max(fittedGLM.linearfilters.Coupling.Filter{pair});
            end
        end
    end
    total_cell_count = cell_count; 
end

%%
default_colors = get(gca,'ColorOrder');
hold on
for cell_count = 1:103
   plot(coupling{cell_count}.CP(:,2),coupling{cell_count}.CP(:,3),'.', 'Color', default_colors(1,:))
end

%%
if 0
for cell_count = 40:142
   plot(coupling{cell_count}.CP(:,2),coupling{cell_count}.CP(:,3),'.', 'Color', default_colors(5,:))
end
for cell_count = 143:210
   plot(coupling{cell_count}.CP(:,2),coupling{cell_count}.CP(:,3),'.', 'Color', default_colors(3,:))
end
for cell_count = 211:265
   plot(coupling{cell_count}.CP(:,2),coupling{cell_count}.CP(:,3),'.', 'Color', default_colors(4,:))
end
end

lim = 3.5;
plot([0 lim],[0 lim])
axis([0 lim 0 lim])
axis square
xlabel('WN Coupling')
ylabel('NSEM Coupling')
title('ON-ON Coupling')

%%

if 0
datarun = load_data('2012-09-27-3/data003');
datarun = load_neurons(datarun);
datarun = load_params(datarun);
plot_rf_fit(datarun, 'On Parasol')
hold on

for i = 40:142
    cid = find(datarun.cell_ids == coupling{i}.info(1));
    start = datarun.vision.sta_fits{cid}.mean;
    for i_pair = 1:6
        cid = find(datarun.cell_ids == coupling{i}.CP(i_pair,1));
        finish = datarun.vision.sta_fits{cid}.mean;
        plot([start(1) finish(1)],[start(2) finish(2)],'k','LineWidth', (0.1+abs(coupling{i}.CP(i_pair,3)))); %-coupling{i}.CP(i_pair,2)
    end
end
end







