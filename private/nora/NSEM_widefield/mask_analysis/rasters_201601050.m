%% Dataruns: run first for everything

clear
regA{1} = 'data003'; cells_regA{1} = [1923;3184;5792;7250];
sigma2{1} = 'data006'; cells_s{1}{1} = [1930;3183;5794;7255];
regB{1} = 'data007'; cells_regB{1} = [1924;3167;5912;7248];
sigma4{1} = 'data008';cells_s{1}{2} = [1922;3181;5914;7246];
cells{1} = 'Maskin2';

regA{2} = 'data007'; cells_regA{2} = [275;2343;3483;4420;5717];
sigma2{2} = 'data010'; cells_s{2}{1} = [345;2464;3379;4413;5716];
regB{2} = 'data011'; cells_regB{2} = [275;2341;3377;4412;5719];
sigma4{2} = 'data012'; cells_s{2}{2} = [202;2465;3481;4411;5717];
cells{2} = 'Maskin';

classification = 'data001';

Analysis_Path = '/Volumes/Analysis/2016-01-05-0/';
fig_save = '/Users/Nora/Desktop/Fig_Output/2016-01-05-0/mVision' ;
mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 15 4])

% classification run
% class_datarun = load_data(['/Volumes/Analysis/2016-01-05-0/streamed/' classification '/' classification]);
% class_datarun = load_params(class_datarun);
% class_datarun = load_neurons(class_datarun);

GS = 0;


%% check stability
%{
for maskin = 2;
    
    datarun = load_data([ Analysis_Path regA{maskin} '-mVision/' regA{maskin} '/' regA{maskin}]);
    datarun = load_neurons(datarun);
    regA_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells_regA{maskin}, 'visual_check', 0);
    datarun = load_data([Analysis_Path regB{maskin} '-mVision/' regB{maskin} '/' regB{maskin}]);
    datarun = load_neurons(datarun);
    regB_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells_regB{maskin}, 'visual_check', 0);
    
    for i = 1:4
        try
            figure;
            IDP_plot_PSTH(regA_data, i);
            hold on; IDP_plot_PSTH(regB_data, i);
        catch
        end
    end
    
end
%}
%% look at all the rasters!

% first 30 are mask
% next 30 are LES
% next 30 are complement
% final 30 are complement + LES
% Load the data
for i = 1:2
    for sigma = 1:2
        if sigma==1
            datarun = load_data([Analysis_Path sigma2{i} '-mVision/' sigma2{i} '/' sigma2{i}]);
            fig_offset = 0;
        else
            datarun = load_data([Analysis_Path sigma4{i} '-mVision/' sigma4{i} '/' sigma4{i}]);
            fig_offset = 3;
        end
        
        datarun = load_neurons(datarun);
        sigma_data = interleaved_data_prep(datarun, 1100, 120, 'cell_spec', cells_s{i}{sigma}, 'visual_check', 0);
        clear datarun
    % separate the conditions
    mask_data{sigma}.testspikes = sigma_data.testspikes(1:30, :);
    %LES_data.testspikes = sigma_data.testspikes(31:60, :);
    comp_data{sigma}.testspikes = sigma_data.testspikes(61:90, :);
    %compLES_data.testspikes = sigma_data.testspikes(91:end, :);
    end
    
    % reg data
    %         datarun = load_data([ Analysis_Path regA{i} '/' regA{i}]);
    %         datarun = load_neurons(datarun);
    %         regA_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells{i}, 'datarun_class', class_datarun, 'visual_check', 0);
    datarun = load_data([ Analysis_Path regB{i} '-mVision/' regB{i} '/' regB{i}]);
    datarun = load_neurons(datarun);
    regB_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells_regB{i}, 'visual_check', 0);
    
    % plot for each cell
    for cell = 1:size(mask_data{1}.testspikes,2)
        figure(2); hold on 
        figure(1); hold on
        reg = IDP_plot_PSTH(regB_data, cell, 'color', 0);
        for sigma=1:2
           color = [1 1 1]*(1-sigma/12)-0.2;
           figure(1); mask = IDP_plot_PSTH(mask_data{sigma}, cell,'color',color);
           %NSEM_Corr(cell, sigma) = err(reg, mask);
           figure(2); comp = IDP_plot_PSTH(comp_data{sigma}, cell,'color', color);
           %surround_struct(cell, sigma) = std(comp);
        end
        figure(1);
        reg = IDP_plot_PSTH(regB_data, cell, 'color', [0 0 0]);
        exportfig(gcf, [fig_save '/maskin' num2str(mod(i,2)+1) '_mask_cell_' num2str(cell)], 'Bounds', 'loose', 'Color', 'rgb');
        clf;
        figure(2);
        exportfig(gcf, [fig_save '/maskin' num2str(mod(i,2)+1) '_comp_cell_' num2str(cell)], 'Bounds', 'loose', 'Color', 'rgb');
        clf;
        
    end
    %save(['SpotStats2016-01-05-0' num2str(i)], 'NSEM_Corr', 'surround_struct');
end

