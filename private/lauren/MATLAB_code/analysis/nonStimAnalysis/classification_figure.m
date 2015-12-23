clear all
% 
% lineColors(1,:) = [90 156 0]/255; %pale grass
% lineColors(2,:) = [255 124 59]/255; %salmon
% lineColors(3,:) = [101 52 255]/255; %purple
% lineColors(4,:) = [52 198 247]/255; %aqua
% lineColors(5,:) = [238 55 128]/255; %calm magenta

blue = [50 70 247]/255;
rust = [.8 .05 0.05];
grass = [90 156 0]/255;


colors(1,:) = 1-(1-blue)*0.7; %ON parasol
colors(2,:) = blue*0.7; %OFF parasol
colors(3,:) = 1-(1-rust)*0.7; %ON midget
colors(4,:) = rust*0.7; %OFF midget
colors(5,:) = grass;


%% datarun for other cell types

datarun = load_data('2010-10-18-3/data009-serial/data009-serial');
datarun.names.obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2010-10-18-3/data009-serial';
goodCells{1} = [10303 10724 20036 20501 50770 60181 70871 90410 100770 130392 180886 270586]; %ON parasol
goodCells{2} = [10107 10543 20212 20586 20663 30423 30691 40036 40649 40769 50426 60451 60766 70006 90902 100226 110272 180136]; %OFF parasol


goodCells{3} = [40061 40633 40650 40933 50139 60151 60211 60886 70331....
    80376 80784 90857 90916 100018 100169 100289 110078 110406 120467....
    120632 130274 130827 160064 160588 170243 170588 190169 190736 190901....
    210196 210421 250001 270076 330139]; %ON midget, excluding really crappy STAs (10621 and 200934)
%goodCells{3} = [10621 40061 40633 40650 40933 50139 60151 60211 60886...
%70331 80376 80784 90857 90916 100018 100169 100289 110078 110406 120467...
%120632 130274 130827 160064 160588 170243 170588 190169 190736 190901...
%200934 210196 210421 250001 270076 330139]; %ON midget, including crappy STAs

goodCells{4} = [10697 30107 40140 50376 50424 50650 60212 60425 70169....
    70529 70739 80396 80787 80797 80934 90170 90230 90303 100827 120093....
    120437 120545 140512 150483 160287 160422 170062 180081 180091 180272....
    180289 200467 210736 230302 240137 240783 280002 300047 300906]; %OFF midget, excluding really crappy STAs (100301)
%goodCells{4} = [10697 30107 40140 50376 50424 50650 60212 60425 70169...
%70529 70739 80396 80787 80797 80934 90170 90230 90303 100301 100827...
%120093 120437 120545 140512 150483 160287 160422 170062 180081 180091...
%180272 180289 200467 210736 230302 240137 240783 280002 300047 300906]; %OFF midget, including crappy STAs

%goodCells{5} = [10873 30573 60724 120197 140921]; %SBCs
    
datarun = load_index(datarun);
datarun = load_neurons(datarun);
%datarun = load_params(datarun, struct('verbose',1,'sync_cell_ids',false));
%datarun = load_sta(datarun,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1,'sync_cell_ids',false);
datarun = load_obvius_sta_fits(datarun, struct('sort_ids', true));
datarun = get_sta_fits_from_obvius(datarun, [goodCells{1} goodCells{2} goodCells{3} goodCells{4}]);

%datarun = get_sta_fits_from_vision(datarun,'all');

% compute the transformation using the clicked points in datarun.piece.array
datarun.stimulus.monitor_x = 640;
datarun.stimulus.monitor_y = 480;
datarun = compute_monitor_to_array_transformation(datarun);

%% dataruns for SBCs

goodCells{5} = [170 588 695 873]; %SBCs ids in non-serial spike sorting
consistentCells = [695 873]; %cells that have qualitatively similar fits in blurred and unblurred cases

%obvius fits to raw STAs
datarun_unblur = load_data('2010-10-18-3/data009-lh/data009-lh');
datarun_unblur = load_index(datarun_unblur);
datarun_unblur = load_neurons(datarun_unblur);
datarun_unblur = load_params(datarun_unblur, struct('verbose',1));
datarun_unblur = load_sta(datarun_unblur,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
datarun_unblur.names.obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2010-10-18-3/data009-lh';
datarun_unblur = load_obvius_sta_fits(datarun_unblur);
datarun_unblur = get_sta_fits_from_obvius(datarun_unblur, consistentCells);

%obvius fits to STAs blurred by Gaussian filter with SD = 0.8 stixels
datarun_blur = load_data('2010-10-18-3/data009-lh-blur-0_8/data009-lh');
datarun_blur = load_index(datarun_blur);
datarun_blur = load_neurons(datarun_blur);
datarun_blur = load_params(datarun_blur, struct('verbose',1));
datarun_blur = load_sta(datarun_blur,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
datarun_blur.names.obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2010-10-18-3/data009-lh-blur-0_8';
datarun_blur = load_obvius_sta_fits(datarun_blur);
datarun_blur = get_sta_fits_from_obvius(datarun_blur, goodCells{5});

%calculate mean area of fits consistent cells for blurred and unblurred cases
original_fits = [datarun_unblur.obvius.sta_fits{get_cell_indices(datarun_unblur, consistentCells)}];
original_fits = [original_fits.sd]; %returns vector of SD pairs: [cell1SD1 cell1SD2 cell2SD1 cell2SD2...]
original_mean_area = pi*mean(original_fits(1:2:end).*original_fits(2:2:end));% area of ellipes = pi*SD1*SD2

blurred_fits = [datarun_blur.obvius.sta_fits{get_cell_indices(datarun_blur, consistentCells)}];
blurred_fits = [blurred_fits.sd];
blurred_mean_area = pi*mean(blurred_fits(1:2:end).*blurred_fits(2:2:end));

%how much to scale blurred SD values by to achieve same mean area of
%consistent cell fits in blurred and unblurred cases
scale_factor = sqrt(original_mean_area/blurred_mean_area);

%tweak blurred fit diameters to have equal area to unblurred version (based
%on consistent cells)
for ii = 1:length(datarun_blur.stas.fits)
    if ~isempty(datarun_blur.stas.fits{ii})
        datarun_blur.obvius.sta_fits{ii}.sd = scale_factor*datarun_blur.obvius.sta_fits{ii}.sd;
        datarun_blur.stas.fits{ii}.sd = scale_factor*datarun_blur.stas.fits{ii}.sd;
    end
end

clear datarun_unblur


% compute the transformation using the clicked points in datarun.piece.array
datarun_blur.stimulus.monitor_x = 640;
datarun_blur.stimulus.monitor_y = 480;
datarun_blur = compute_monitor_to_array_transformation(datarun_blur);


% %% separate out cells that lie over the array
% 
% %midgets and parasols
% for ii = 1:4
%     for jj = 1:length(goodCells{ii})
%         cellInfo(jj).id = goodCells{ii}(jj);
%         cellInfo(jj).type = 'onPar'; %not true for all ii but doesn't matter for this purpose
%     end
%     [goodCellsOverArray{ii} goodCellsOffArray{ii}] = removeSmallSigEdgeCells(cellInfo, datarun.names.rrs_ei_path, 'cutOffMult', 0.5);
%     clear cellInfo
% end
% 
% 
% %SBCs
% ii = 5;
% for jj = 1:length(goodCells{ii})
%     cellInfo(jj).id = goodCells{ii}(jj);
%     cellInfo(jj).type = 'onPar'; %not true for all ii but doesn't matter for this purpose
% end
% [goodCellsOverArray{ii} goodCellsOffArray{ii}] = removeSmallSigEdgeCells(cellInfo, datarun_blur.names.rrs_ei_path, 'cutOffMult', 0.5);
% clear cellInfo
% 
% for ii = 1:5
%     goodCells{ii} = [goodCellsOverArray{ii}.id];
% end
    
%% plot full mosaics

figure
plot_rf_summaries(datarun, goodCells{1}, 'fit_width', 1, 'array', true, 'fit_color', colors(1,:))
hold on
plot_rf_summaries(datarun, goodCells{2}, 'fit_width', 1, 'fit_color', colors(2,:))
plot_rf_summaries(datarun, goodCells{3}, 'fit_width', 1, 'fit_color', colors(3,:))
plot_rf_summaries(datarun, goodCells{4}, 'fit_width', 1, 'fit_color', colors(4,:))
%plot_rf_summaries(datarun, goodCells{5}, 'fit_width', 1, 'fit_color', lineColors(5,:))
plot_rf_summaries(datarun_blur, goodCells{5}, 'fit_width', 1, 'fit_color', colors(5,:))



%% get sta timecourses for each cell

included_cells = [];

    
%midgets and parasols
paramsFile = edu.ucsc.neurobiology.vision.io.ParametersFile(datarun.names.rrs_params_path);
for ii=1:length(datarun.cell_ids)
    
    t = paramsFile.getArrayCell(datarun.cell_ids(ii), 'RedTimeCourse');
    t1 = paramsFile.getArrayCell(datarun.cell_ids(ii), 'GreenTimeCourse');
    t2 = paramsFile.getArrayCell(datarun.cell_ids(ii), 'BlueTimeCourse');
    %datarun.green_time_course{i} = t/sum(abs(t));
    %datarun.green_time_course{i} = t/sum(abs(t+t1+t2));
    
    rms = sqrt(sum(t.^2) + sum(t1.^2) + sum(t2.^2)); %normalization: rms of all 3 STA primaries
    
    datarun.red_time_course{ii} = t/rms; 
    datarun.green_time_course{ii} = t1/rms;
    datarun.blue_time_course{ii} = t2/rms; 
    
    
    %datarun.green_time_course{i} = t/sqrt(sum((t+t1+t2).^2));
    %t = paramsFile.getArrayCell(datarun.cell_ids(ii), 'Auto');
    %datarun.auto{i} = t/sum(abs(t));
    %datarun.auto{ii} = t/sqrt(sum(t.^2)); %normalizes acf by rms
    %datarun.auto{i} = t;
end
paramsFile(close)


%SBCs
paramsFile_blur = edu.ucsc.neurobiology.vision.io.ParametersFile(datarun_blur.names.rrs_params_path);
for ii=1:length(datarun_blur.cell_ids)
    t = paramsFile_blur.getArrayCell(datarun_blur.cell_ids(ii), 'RedTimeCourse');
    t1 = paramsFile_blur.getArrayCell(datarun_blur.cell_ids(ii), 'GreenTimeCourse');
    t2 = paramsFile_blur.getArrayCell(datarun_blur.cell_ids(ii), 'BlueTimeCourse');
    
    rms = sqrt(sum(t.^2) + sum(t1.^2) + sum(t2.^2)); %normalization: rms of all 3 STA primaries
    
    datarun_blur.red_time_course{ii} = t/rms; 
    datarun_blur.green_time_course{ii} = t1/rms;
    datarun_blur.blue_time_course{ii} = t2/rms; 
end
paramsFile_blur(close)


%% plot STA timecourses and collect in cell array sta_all

sta_all = cell(1,5);

%midgets and parasols
for ii = 1:4
    indices = get_cell_indices(datarun, goodCells{ii});
    sta_all{ii} = [];
    figure('position', [100 100 120 200]); hold on
    for jj = 1:length(indices)
        plot(datarun.red_time_course{indices(jj)}(end-12:end),'r-')
        plot(datarun.green_time_course{indices(jj)}(end-12:end),'g-')
        plot(datarun.blue_time_course{indices(jj)}(end-12:end),'b-')
        
        sta_all{ii} = cat(3, sta_all{ii}, [datarun.red_time_course{indices(jj)} datarun.green_time_course{indices(jj)} datarun.blue_time_course{indices(jj)}]);
    end
    sta_mean = mean(sta_all{ii}, 3);
    
    plot(sta_mean(end-12:end,1), 'color', [0.8 0 0], 'LineWidth', 2)
    plot(sta_mean(end-12:end,2), 'color', [0 0.8 0], 'LineWidth', 2)
    plot(sta_mean(end-12:end,3), 'color', [0 0 0.8], 'LineWidth', 2)
    set(gca, 'xlim', [1 13], 'ylim', [-1 1], 'xtick', [1 7 13], 'ytick', [], 'xticklabel', [-200 -100 0])    
    
    hold off
end

%SBCs
ii = 5;
indices = get_cell_indices(datarun_blur, goodCells{ii});
sta_all{ii} = [];
figure('position', [100 100 120 200]); hold on
for jj = 1:length(indices)
    plot(datarun_blur.red_time_course{indices(jj)}(end-12:end),'r-')
    plot(datarun_blur.green_time_course{indices(jj)}(end-12:end),'g-')
    plot(datarun_blur.blue_time_course{indices(jj)}(end-12:end),'b-')
    
    sta_all{ii} = cat(3, sta_all{ii}, [datarun_blur.red_time_course{indices(jj)} datarun_blur.green_time_course{indices(jj)} datarun_blur.blue_time_course{indices(jj)}]);
end
sta_mean = mean(sta_all{ii}, 3);

plot(sta_mean(end-12:end,1), 'color', [0.8 0 0], 'LineWidth', 2)
plot(sta_mean(end-12:end,2), 'color', [0 0.8 0], 'LineWidth', 2)
plot(sta_mean(end-12:end,3), 'color', [0 0 0.8], 'LineWidth', 2)
set(gca, 'xlim', [1 13], 'ylim', [-1 1], 'xtick', [1 7 13], 'ytick', [], 'xticklabel', [-200 -100 0])

hold off


%%
% make plot based on visual parameters

%calculate STA PC and rf diameters
c=[];
id=[];
sta_data = [];
rf_diam = [];
sta_in_class = [];

%midgets and parasols
for ii = 1:4
    indices = get_cell_indices(datarun, goodCells{ii});
    for jj=1:length(indices)
        id=[id datarun.cell_ids(indices(jj))];
        c=[c; colors(ii,:)]; %stores color of each point
        
        sta_data = [sta_data; [datarun.green_time_course{indices(jj)}' datarun.red_time_course{indices(jj)}' datarun.blue_time_course{indices(jj)}']];
        rf_diam = [rf_diam sqrt(datarun.stas.fits{indices(jj)}.sd(1)*datarun.stas.fits{indices(jj)}.sd(2))]; %equivalent circular radius, relative
        
        sta_in_class = [sta_in_class; [datarun.red_time_course{indices(jj)}' datarun.green_time_course{indices(jj)}' datarun.blue_time_course{indices(jj)}']];
        
        %r=[r ext(datarun.green_time_course{indices(jj)})]; %ext returns most extreme value and index of most extreme value
        %rr=[rr pi*prod(datarun.vision.sta_fits{indices(jj)}.sd)*8*5.75];
        
    end
end

%SBCs
ii=5;
indices = get_cell_indices(datarun_blur, goodCells{ii});
for jj=1:length(indices)
    id=[id datarun_blur.cell_ids(indices(jj))];
    c=[c; colors(ii,:)]; %stores color of each point
    
    sta_data = [sta_data; [datarun_blur.green_time_course{indices(jj)}' datarun_blur.red_time_course{indices(jj)}' datarun_blur.blue_time_course{indices(jj)}']];
    rf_diam = [rf_diam sqrt(datarun_blur.stas.fits{indices(jj)}.sd(1)*datarun_blur.stas.fits{indices(jj)}.sd(2))]; %equivalent circular radius, relative
    
    sta_in_class = [sta_in_class; [datarun_blur.red_time_course{indices(jj)}' datarun_blur.green_time_course{indices(jj)}' datarun_blur.blue_time_course{indices(jj)}']];
    
    %r=[r ext(datarun_blur.green_time_course{indices(jj)})]; %ext returns most extreme value and index of most extreme value
    %rr=[rr pi*prod(datarun_blur.vision.sta_fits{indices(jj)}.sd)*8*5.75];
    
end


[coeff pc] = princomp(sta_data);

pc1 = pc(:,1);

%%

figure('position', [100 100 300 400])

hold on
for ii = 1:length(pc1)
    plot(pc1(ii), rf_diam(ii), '.', 'MarkerFaceColor', c(ii,:), 'MarkerEdgeColor', c(ii,:), 'MarkerSize', 12)
end
set(gca,'PlotBoxAspectRatio',[0.6 1 1],'Xtick',[0],'Ytick',[0]);
xlabel('STA PC1')
ylabel('RF diameter')
set(gca, 'xlim', [min(pc1) max(pc1)] + [-0.1 0.1]*range(pc1), 'ylim', [min(rf_diam) max(rf_diam)] + [-0.1 0.1]*range(rf_diam))
box on

