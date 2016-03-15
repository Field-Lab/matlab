%% Dataruns
clear
regA = 'data005';
mask = 'data007';
LES = 'data008';
comp = 'data009';
compLES = 'data010';
% regB = 'data010';
cells = 'Maskin';

classification = 'data002';

Analysis_Path = '/Volumes/Analysis/2015-12-18-2/map-from-data002-stream/';
fig_save = '/Users/Nora/Desktop/Fig_Output';
mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')

% classification run
class_datarun = load_data(['/Volumes/Analysis/2015-12-18-2/streamed/' classification '/' classification]);
class_datarun = load_params(class_datarun);
class_datarun = load_neurons(class_datarun);

%% check stability

datarun = load_data([ Analysis_Path regA '/' regA]);
datarun = load_neurons(datarun);
regA_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);
datarun = load_data([Analysis_Path regB '/' regB]);
datarun = load_neurons(datarun);
regB_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);

for i = 1:6
    try
        figure;
        IDP_plot_PSTH(regA_data, i);
        hold on; IDP_plot_PSTH(regB_data, i);
    catch
    end
end

%% load up the conditions
datarun = load_data([ Analysis_Path mask '/' mask]);
datarun = load_neurons(datarun);
mask_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);

datarun = load_data([ Analysis_Path LES '/' LES]);
datarun = load_neurons(datarun);
LES_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);

datarun = load_data([ Analysis_Path comp '/' comp]);
datarun = load_neurons(datarun);
comp_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);

datarun = load_data([ Analysis_Path compLES '/' compLES]);
datarun = load_neurons(datarun);
compLES_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);

datarun = load_data([ Analysis_Path regB '/' regB]);
datarun = load_neurons(datarun);
regA_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);

%% plot for each cell
for cell = 1:size(regA_data.testspikes,2)
    
    save_name = [fig_save '/sigma' num2str(1) '_cell' num2str(cell) '_'];
    
    % reg versus mask
    figure; hold on;
    reg = IDP_plot_PSTH(regA_data, cell);
    mask = IDP_plot_PSTH(mask_data, cell);
    legend('Full', 'Masked')
    print('-dpdf', [save_name 'reg_mask'])
    close(gcf)
    
    % reg versus complement
    figure; hold on;
    plot(reg)
    comp = IDP_plot_PSTH(comp_data, cell);
    legend('Full', 'Complement')
    print('-dpdf', [save_name 'reg_comp'])
    close(gcf)
    
    % reg versus complement + mask
    figure; hold on;
    plot(reg)
    plot(mask+comp)
    set(gcf, 'Position', [ 500 500 800 200])
    legend('Full', 'Masked + Comp')
    print('-dpdf', [save_name 'reg_mask_comp'])
    close(gcf)
    
    % mask versus LES
    figure; hold on;
    plot(mask)
    LES = IDP_plot_PSTH(LES_data, cell);
    legend('Masked', 'LES')
    print('-dpdf', [save_name 'mask_LES'])
    close(gcf)
    
    % reg versus LES+complement
    figure; hold on;
    plot(reg)
    plot(LES+comp)
    set(gcf, 'Position', [ 500 500 800 200])
    compLES = IDP_plot_PSTH(compLES_data, cell);
    legend('Full', 'LES + Complement', 'LES and Complement (run together)')
    print('-dpdf', [save_name 'full_compLES'])
    close(gcf)
end

