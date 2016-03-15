%% dataset differences quantified!!!

clear
tic
% 1-3 exp, 1-2 On OFF
% default_colors = get(gca,'ColorOrder');
% ExpColors = [default_colors(2,:); default_colors(5,:); default_colors(1,:)];
% TypeSymbols = {'+','.'};
% TypeLines = {'--', ''};

datasets = {'2012-08-09-3', '2012-09-27-3','2013-08-19-6'};
cell_types = {'On Parasol', 'Off Parasol'};

% fractional variance
base_path = '/Volumes/Lab/Users/akheitman/NSEM_Home/Cell_Selection/Raster_Metrics/base_crm_findPS_';
hold on
for data = 1:3
    load([base_path datasets{data} '.mat'])
    for cell_type = 1:2
        WN_scores{data, cell_type} = raster_scores.celltype{cell_type}.scores_WN.glm_bps./raster_scores.celltype{cell_type}.scores_WN.crm_bps;
        WN_scores{data, cell_type}(WN_scores{data, cell_type}<0) = 0;
        NSEM_scores{data, cell_type} = raster_scores.celltype{cell_type}.scores_NSEM.glm_bps./raster_scores.celltype{cell_type}.scores_NSEM.crm_bps;
        NSEM_scores{data, cell_type}(NSEM_scores{data, cell_type}<0) = 0;
        if cell_type == 1
            cell_ids{data, cell_type} = raster_scores.ONP;
        else
            cell_ids{data, cell_type} = raster_scores.OFFP;
        end
    end
end
figure; 
exp3plot(WN_scores,NSEM_scores)
exp3_legend(gca, 1)
xlabel('WN BPS')
ylabel('NSEM BPS')
xlim([0 1])
ylim([0 1])
axis square
plot([0 1], [0 1], 'k')
clear raster_scores i cell_type base_path data
figure; 
exp3plot(WN_scores, 1);
ylim([0 1])
exp3plot(NSEM_scores, 1);
ylim([0 1])

%% VISION STATS AND GENERATOR SIGNAL

for data = 1:length(datasets)
    
    exp_info = experimentinfoNSEM(datasets{data});
    datarun = load_data([datasets{data} '/' exp_info.dr.mas]);
    datarun = load_params(datarun);
    datarun = load_neurons(datarun);
    
    for cell_type = 1:length(cell_types)
        cells = get_cell_indices(datarun, cell_ids{data, cell_type});
        
        % average timecourse and RGB
        timecourse{data,cell_type} = mean([datarun.vision.timecourses(cells).r datarun.vision.timecourses(cells).g datarun.vision.timecourses(cells).b]');
        if data == 1
           temp = timecourse{data,cell_type}(1:2:end);
           timecourse{data,cell_type}(1:15) = zeros(15,1);
           timecourse{data,cell_type}(16:end) = temp;
           clear temp
        end
        %timecourse{data,cell_type} = s(:,1);
        RGB_avg = zeros(3,1);
        for cell = 1:length(cells)
            RGB_avg = RGB_avg + abs(RGB_weights(datarun, cells(cell)));
            sta_sd = datarun.vision.sta_fits{cells(cell)}.sd;
            %sta_ell{data,cell_type}(cell) = abs(diff(sta_sd));
            sta_mean{data,cell_type}(cell) = mean(sta_sd);
            spike_rate{data,cell_type}(cell) = length(datarun.spikes{cells(cell)})/datarun.duration;
        end
        RGB_avg = RGB_avg/length(cells);
        RGB_avg = RGB_avg/max(RGB_avg);
        RGB{data,cell_type} = RGB_avg; 
    end
end

sta_mean{1,1} = sta_mean{1,1}*0.8;
sta_mean{1,2} = sta_mean{1,2}*0.8;

exp3plot(timecourse, 0) % need to fix this
exp3plot(RGB, 0)
exp3plot(sta_mean,1)
exp3plot(spike_rate,1)
exp3plot(NSEM_scores, spike_rate)
clear timecourse RGB sta_mean spike_rate sta_ell sta_sd RGB_avg exp_info datarun data cell_type cells

%% GLM FITS

%datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_noCP_p8IDp8_shortlist/standardparams/';
datapath='/Volumes/Lab/Users/akheitman/NSEM_Home/GLM_Output_Analysis/fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams/';

% exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{1}='WN_mapPRJ/';
fittypepath{2}='NSEM_mapPRJ/';

cell_type_str = {'ONPar_', 'OFFPar_'};

for data=1:3
    for cell_type = 1:2
        n_cells = length(cell_ids{data, cell_type});
        
        %Initialize matrices
        NSEM_mu{data,cell_type}=zeros(n_cells,1);
        WN_time_2_peak{data,cell_type}=zeros(n_cells,1);
        NSEM_PS{data,cell_type}=zeros(n_cells,120);
        NSEM_K_space{data,cell_type}=zeros(n_cells, 13, 13);
        NSEM_K_time{data,cell_type}=zeros(n_cells, 30);
        WN_mu{data,cell_type}=zeros(n_cells,1);
        WN_PS{data,cell_type}=zeros(n_cells,120);
        WN_K_space{data,cell_type}=zeros(n_cells, 13, 13);
        WN_K_time{data,cell_type}=zeros(n_cells, 30);
        
        for cell=1:n_cells
            %try
            load([datapath fittypepath{2} datasets{data} '/' cell_type_str{cell_type} num2str(cell_ids{data,cell_type}(cell)) '.mat']);
            NSEM_K_space{data,cell_type}(cell,:,:)=fittedGLM.linearfilters.Stimulus.space_rk1;
            NSEM_K_time{data,cell_type}(cell,:)=fittedGLM.linearfilters.Stimulus.time_rk1;
            NSEM_mu{data,cell_type}(cell)=fittedGLM.linearfilters.TonicDrive.Filter;
            NSEM_PS{data,cell_type}(cell,:)=fittedGLM.linearfilters.PostSpike.Filter;
            load([datapath fittypepath{1} datasets{data} '/' cell_type_str{cell_type} num2str(cell_ids{data,cell_type}(cell)) '.mat']);
            WN_K_space{data,cell_type}(cell,:,:)=fittedGLM.linearfilters.Stimulus.space_rk1;
            WN_K_time{data,cell_type}(cell,:)=fittedGLM.linearfilters.Stimulus.time_rk1;
            WN_mu{data,cell_type}(cell)=fittedGLM.linearfilters.TonicDrive.Filter;
            WN_PS{data,cell_type}(cell,:)=fittedGLM.linearfilters.PostSpike.Filter;
            [~,WN_time_2_peak{data,cell_type}(cell)] = max(abs(WN_K_time{data,cell_type}(cell,:)));
            
            %catch
            %    disp(['No fit for cell' num2str(cell_ids{data,cell_type}(cell))])
            %end
        end
    end
end

figure; exp3plot(NSEM_K_time,0)
figure; exp3plot(WN_K_time,0)
% figure; exp3plot(NSEM_PS,0)
% figure; exp3plot(WN_PS,0)
% figure; exp3plot(WN_scores, WN_mu)
% figure; exp3plot(NSEM_scores, NSEM_mu)
% figure; exp3plot(WN_scores, WN_time_2_peak)

clear n_cells fittypepath data datapath cell cell_type cell_type_str fittedGLM NSEM_K_space NSEM_K_time NSEM_mu NSEM_PS WN_K_space WN_K_time WN_mu WN_PS WN_time_2_peak


%% SPIKE SORTING!

tic
for fit=1
    for data = 1:3
        
        exp_info = experimentinfoNSEM(datasets{data});
        
        if fit == 1
            datarun = [exp_info.dr.slvWN '-from-' exp_info.dr.mas];
        else
            datarun = [exp_info.dr.slvNSEM '-from-' exp_info.dr.mas];
        end
        stats_file = ['/Volumes/Analysis/' datasets{data} '/' datarun '/' datarun '-stats.txt'];
        %model_file = ['/Volumes/Analysis/' datasets{data} '/' datarun '/' datarun '.model'];
        stats = load(stats_file);
        cell_id = stats(:,4);
        
        for cell_type = 1:2
            disp([data cell_type])
            toc
            for cell = 1:length(cell_ids{data, cell_type})
                
                %index = find(cell_id == cell_ids{data, cell_type}(cell));
                all_stats{data, cell_type}(cell,:) = stats(cell_id == cell_ids{data, cell_type}(cell),:);
                %contamination_index{data, cell_type}(cell) = stats(index,9);
                %best_union{data, cell_type}(cell) = stats(index,13);
                %worst_likelihood{data, cell_type}(cell) = stats(index,14);
                
                %electrode = stats(index,5);
                %cluster_number = stats(index,6);
                %model = edu.ucsc.neurobiology.vision.io.ClusteringModelFile(model_file);
                % get the clusters on electrode 1
                %clusters = model.getNeuronExtraction(electrode);
                %cluster_prob{data, cell_type}(cell) = clusters.probability(cluster_number+1);
                %model.close
                %clear clusters
            end
        end
        
    end
%     if fit == 1
%         WN_CI = contamination_index;
%         WN_BU = best_union;
%         WN_WL = worst_likelihood;
%         WN_CP = cluster_prob;
%     else
%         NSEM_CI = contamination_index;
%         NSEM_BU = best_union;
%         NSEM_WL = worst_likelihood;
%         NSEM_CP = cluster_prob;
%     end
    clear contamination_index best_union worst_likelihood cluster_prob
end
toc

%%
for s = 10
    for i = 1:3
        for j = 1:2
            temp1{i,j} = all_stats{i,j}(:,s);
            %temp2{i,j} = NSEM_stats{i,j}(:,s);
        end
    end
    figure(s);
   exp3plot(temp1, WN_scores)
   %axis square
   %hold on; plot([0 10], [0 10], 'k'); hold off;
end

%%
default_colors = get(gca,'ColorOrder');
ExpColors = [default_colors(2,:); default_colors(5,:); default_colors(1,:)];
hold on
count = 1;
hist_counts = zeros(6,2);
for s = 11
    for i = 1:3
        for j = 1:2
            %subplot(3,2,count)
            hist_counts(count,:) = histcounts(all_stats{i,j}(:,s), 2, 'Normalization', 'probability');
            %hist_counts = histogram(all_stats{i,j}(:,s), 2, 'Normalization', 'probability', 'FaceColor', ExpColors(i,:));
            %temp1{i,j} = mean(all_stats{i,j}(:,s));
            %temp2{i,j} = mean(WN_scores{i,j});
            %temp2{i,j} = NSEM_stats{i,j}(:,s);
            %plot(j+2*i-1.5, temp1{i,j}, 'o')
                
            count = count +1;
            %set(gca, 'XTick', [0.25 0.75]);
            %set(gca, 'XTickLabels', {'No', 'Yes'})
            %ylim([0 1])
        end
    end

   figure(s);
   
   x = {'On', 'Off', 'On', 'Off', 'On', 'Off'};
   b = bar(hist_counts(:,2));
   disp('here')
   set(b(1), 'FaceColor', 'none');
   set(b(2), 'FaceColor', default_colors(2,:));
   set(b(3), 'FaceColor', ExpColors(2,:));
   set(b(4), 'FaceColor', ExpColors(2,:));
   set(b(5), 'FaceColor', ExpColors(3,:));
   set(b(6), 'FaceColor', ExpColors(3,:));
   set(gca, 'XTickLabels', x)
   
   %axis square
   %hold on; plot([0 10], [0 10], 'k'); hold off;
end


%%
exp3plot(WN_BU, 1)
exp3plot(WN_CP, 1)
exp3plot(WN_WL, 1)
exp3plot(NSEM_WL, 1)
exp3plot(NSEM_CP, 1)
exp3plot(NSEM_BU, 1)

exp3plot(WN_BU, NSEM_BU)
axis square
plot([0 1], [0 1], 'k')
axis([0 0.5 0 0.5])
ylabel('WN')
xlabel('NSEM')
exp3plot(WN_WL, NSEM_WL)
axis([0 120 0 120])
axis square
plot([0 120], [0 120], 'k')
ylabel('WN')
xlabel('NSEM')
exp3plot(WN_CP, NSEM_CP)
axis([0 1 0 1])
axis square
plot([0 1], [0 1], 'k')
ylabel('WN')
xlabel('NSEM')

%%
exp3plot(WN_BU, WN_scores)
exp3plot(WN_WL, WN_scores)
exp3plot(WN_CP, WN_scores)

clear worst_likelihood datarun best_union cell_id contamination_index index clusters model cluster_number electrode exp_info stats_file stats data fit cell_type cell model_file
clear NSEM_BU NSEM_CI NSEM_CP NSEM_WL WN_BU WN_CI WN_CP WN_WL cluster_prob
%% VIZ
% generator signal (already calculated and saved averages)
hold on
default_colors = get(gca,'ColorOrder');
load('/Users/Nora/Desktop/Generator Signal/avg_GS_2012-08-09-3_On Parasol.mat');
plot(avg_GS(1,:), avg_GS(2,:), '--', 'Color', default_colors(2,:))
load('/Users/Nora/Desktop/Generator Signal/avg_GS_2012-08-09-3_Off Parasol.mat');
plot(avg_GS(1,:), avg_GS(2,:), 'Color', default_colors(2,:), 'LineWidth', 2)

load('/Users/Nora/Desktop/Generator Signal/avg_GS_2012-09-27-3_On Parasol.mat');
plot(avg_GS(1,:), avg_GS(2,:), '--', 'Color', default_colors(5,:))
load('/Users/Nora/Desktop/Generator Signal/avg_GS_2012-09-27-3_Off Parasol.mat');
plot(avg_GS(1,:), avg_GS(2,:), 'Color', default_colors(5,:), 'LineWidth', 2)

load('/Users/Nora/Desktop/Generator Signal/avg_GS_2013-08-19-6_On Parasol.mat');
plot(avg_GS(1,:), avg_GS(2,:), '--', 'Color', default_colors(1,:))
load('/Users/Nora/Desktop/Generator Signal/avg_GS_2013-08-19-6_Off Parasol.mat');
plot(avg_GS(1,:), avg_GS(2,:), 'Color', default_colors(1,:), 'LineWidth', 2)

exp3_legend(gca,1)
title('generator signal')
clear avg_GS


