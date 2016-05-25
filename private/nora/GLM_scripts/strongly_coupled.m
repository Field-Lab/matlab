clear

datapath='/Volumes/Lab/Users/akheitman/NSEM_Home/GLM_Output_Analysis/fixedSP_rk1_linear_MU_PS_CP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{2}='NSEM_mapPRJ/';
fittypepath{1}='WN_mapPRJ/';

cell_count = 0;
total_cell_count = 0;

% coupling.info = cellid, exp, BPS WN, BPS NSEM
% coupling.CP = paired cell, CP_WN, CP_NSEM

for exp = 1:3
    for fittype = 2
        cell_count = total_cell_count;
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) '*.mat']);
        n_cells=length(matfiles);
        
        disp(['Experiment ' num2str(exp) ' Stimulus ' fittypepath{fittype} ' with ' num2str(n_cells) ' Cells'])

        % Load cells
        for file=1:n_cells
            cell_count = cell_count + 1;
            load([datapath fittypepath{1} exp_names(exp,:) matfiles(file).name]);
            coupling{cell_count}.info(1) = str2double(matfiles(file).name(7:(end-4)));
            coupling{cell_count}.CP(:,1) = double(fittedGLM.cellinfo.pairs);
            
            % WN
            coupling{cell_count}.info(2) = fittedGLM.xvalperformance.glm_normedbits;
            for pair = 1:6
                coupling{cell_count}.CP(pair,2) = max(fittedGLM.linearfilters.Coupling.Filter{pair});
            end
            
            % NSEM
            load([datapath fittypepath{2} exp_names(exp,:) matfiles(file).name]);
            coupling{cell_count}.info(3) = fittedGLM.xvalperformance.glm_normedbits;
            for pair = 1:6
                coupling{cell_count}.CP(pair,3) = max(fittedGLM.linearfilters.Coupling.Filter{pair});
            end
            

        end
    end
    total_cell_count = cell_count
end

%%
default_colors = get(gca,'ColorOrder');
hold on
for cell_count = 1:10
   plot(coupling{cell_count}.CP(:,2),coupling{cell_count}.CP(:,3),'.', 'Color', default_colors(1,:))
end

%%
hold on
for cell_count = 1:9
   plot(coupling{cell_count}.CP(:,2),coupling{cell_count}.CP(:,3),'.', 'Color', default_colors(1,:))
end
for cell_count = 10:18
   plot(coupling{cell_count}.CP(:,2),coupling{cell_count}.CP(:,3),'.', 'Color', default_colors(2,:))
end
for cell_count = 19:28
   plot(coupling{cell_count}.CP(:,2),coupling{cell_count}.CP(:,3),'.', 'Color', default_colors(5,:))
end

lim = 2;
plot([0 lim],[0 lim])
axis([0 lim 0 lim])
axis square
xlabel('WN Coupling')
ylabel('NSEM Coupling')
title('ON-ON Coupling')

%%

% datarun = load_data('2012-09-27-3/data003');
% datarun = load_neurons(datarun);
% datarun = load_params(datarun);
% plot_rf_fit(datarun, 'On Parasol')
% hold on
% plot_rf_fit(datarun, [91 1909 2360 6858], 'edge', false, 'fill', true, 'alpha', 0.1)
% 
% for i = 1:length(coupling)
%     cid = find(datarun.cell_ids == coupling{i}.info(1));
%     start = datarun.vision.sta_fits{cid}.mean;
%     for i_pair = 1:6
%         cid = find(datarun.cell_ids == coupling{i}.CP(i_pair,1));
%         finish = datarun.vision.sta_fits{cid}.mean;
%         plot([start(1) finish(1)],[start(2) finish(2)],'k','LineWidth', (0.1+abs(coupling{i}.CP(i_pair,3)))); %-coupling{i}.CP(i_pair,2)
%     end
% end







