clear

piece = '2016-02-17-8';
reg{1} = 'data003';
reg{2} = 'data005';

% Map Types
% 1 = map to stream data
Analysis_Path{1} = ['/Volumes/Analysis/' piece '/map-from-data000/'];
% 2 = concatenated mapping
Analysis_Path{2} = ['/Volumes/Analysis/' piece '/data000-data006/'];
% 3 = separate mVision
Analysis_Path{3} = ['/Volumes/Analysis/' piece '/mVision/'];
% 4 = separate vision
Analysis_Path{4} = ['/Volumes/Analysis/' piece '/'];

% Maskin

%{
masking = 'data002';
%{
% n_masks = 7;
% sigmas = [2 2 4 4 5 5 6 6 4 4 5 5 6 6];
% % subgroup 1
% cell_idx{1} = 1:3;
% mask_conditions{1} = [1 3 5 7];
% comp_conditions{1} = [2 4 6 8];
% % subgroup 1
% cell_idx{2} = 4:6;
% mask_conditions{2} = [1 9 11 13];
% comp_conditions{2} = [2 10 12 14];
% Different cell ids for the four map types
%}
cells{1}{1} = [917, 2731, 3452, 4876, 5973, 7250]; % these will have the same ids
cells{1}{2} = cells{1}{1};
cells{2}{1} = [918;2731;3456;4921;6065;7128]; % cell ids for concatenated mapping, matched to streamed
cells{2}{2} = cells{2}{1};
cells{3}{1} = [916;2732;3453;4923;6062;7141];% cell ids for mVision data003, matched to streamed
cells{3}{2} = [916;2733;3453;4876;5957;7144];% cell ids for mVision data005, matched to streamed
cells{4}{1} = [917;2731;3454;4877;6065;7128];% cell ids for Vision data003, matched to streamed
cells{4}{2} = [917;2733;3453;4923;5958;7129];% cell ids for Vision data005, matched to streamed
%}

piece = '2016-02-17-8';

class = 'data006';
masking = 'data008';
split = 'data001';
reg{1} = 'data010';

% Midgets
cells{1}{1} = [1157, 1561, 2388, 3031, 3977, 4867, 5451, 5956, 6458, 7265];
cells{1}{2} = cells{1}{2};
cells{3}{1} = [916;2732;3453;4923;6062;7141];% cell ids for mVision data005, matched to streamed
cells{3}{2} = [916;2733;3453;4876;5957;7144];% cell ids for mVision data010, matched to streamed
cells{4}{1} = [[1157;1562;1756;766;3977;4746;5508;5956;6451;7263]];% cell ids for Vision data005, matched to streamed
cells{4}{2} = [917;2733;3453;4923;5958;7129];% cell ids for Vision data010, matched to streamed
n_masks{1} = 2;
sigmas{1} = [2 2 4 4];
% subgroup 1
cell_idx{1}{1} = 1:10;
mask_conditions{1}{1} = [1 3];
comp_conditions{1}{1} = [2 4];


fig_save = ['/Users/Nora/Desktop/Fig_Output/' piece];
mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 15 4])

%%

for i_map = 1:4
    
    % load reg repeats
    for i_reg = 1:2
        if i_map == 1
            datarun{i_reg} = load_data([ Analysis_Path{i_map} reg{i_reg} '-cf/' reg{i_reg} '-cf']);
        else
            datarun{i_reg} = load_data([ Analysis_Path{i_map} reg{i_reg} '/' reg{i_reg}]);
        end
        datarun{i_reg} = load_neurons(datarun{i_reg});
        reg_data{i_reg} = interleaved_data_prep(datarun{i_reg}, 1100, 30, 'cell_spec', cells{i_map}{i_reg},'visual_check', 0);
        clear datarun
    end
    
    %plot reg repeats for stability check
    for i_cell = 6%1:size(reg_data{1}.testspikes,2)
        figure; hold on
        for i_reg = 1:2
            reg_PSTH{i_map}{i_reg}{i_cell} = IDP_plot_PSTH(reg_data{i_reg}, i_cell, i_reg);
        end
        title(['cell ' num2str(i_cell) ', map' num2str(i_map)]);
        exportfig(gcf, [fig_save '/' 'cell ' num2str(i_cell) ', map_real_' num2str(i_map)], 'Bounds', 'loose', 'Color', 'rgb')
    end
end

%%

for i_cell = 1:size(reg_data{1}.testspikes,2)
    for i_reg = 1:2
        figure; hold on
        for i_map = 1:4
            plot(reg_PSTH{i_map}{i_reg}{i_cell})
        end
        legend('Mapped from stream, hand clustered', 'Concatenated Mapping', 'EI Mapped, MVision', 'EI Mapped, Vision', 'Location', 'eastoutside')
        title(['cell ' num2str(i_cell) ', map' num2str(i_reg)]);
        xlim([0 500])
        exportfig(gcf, [fig_save '/' 'cell ' num2str(i_cell) ', map' num2str(i_reg)], 'Bounds', 'loose', 'Color', 'rgb')
    end
end

%{
    % replace with cellfindered data if an option
datarun_cf = load_data([ Analysis_Path masking{i_group} '-cf/' masking{i_group} '-cf']);
datarun_cf = load_neurons(datarun_cf);
mask_data_cf = interleaved_data_prep(datarun_cf, 1100, 480, 'cell_spec', 'all', 'visual_check', 0);

for i_cell = 1:length(cells{i_group})
    try
        mask_data.testspikes(:,i_cell) = mask_data_cf.testspikes(:,cells{i_group}(i_cell) == datarun_cf.cell_ids);
    catch
        disp(['cell ' num2str(cells{i_group}(i_cell)) ' not cell findered'])
    end
end
clear datarun_cf mask_data_cf
%}
    

