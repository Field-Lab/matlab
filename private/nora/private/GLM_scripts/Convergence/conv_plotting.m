clear
datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_CP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{2}='NSEM_mapPRJ/';
fittypepath{1}='WN_mapPRJ/';

% Conv_Blocks = [1 2 3 5 7 9 11 20 30 40 50 57];
Conv_Blocks = [11 20 30];

experiment = 1:3;

for fittype=2
    for exp = experiment;
        BPS_temp = zeros(12, length(Conv_Blocks));
        for conv = 1:length(Conv_Blocks)
            
            % Get file list
            matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) 'conv_blocks_' num2str(Conv_Blocks(conv)) '/*.mat']);
            n_cells = length(matfiles);
            
            % Collect info from files
            for file=1:n_cells
                load([datapath fittypepath{fittype} exp_names(exp,:) 'conv_blocks_' num2str(Conv_Blocks(conv)) '/' matfiles(file).name]);
                BPS_temp(file,conv) = fittedGLM.xvalperformance.glm_normedbits;
            end
        end
        BPS{exp,fittype} = BPS_temp(1:n_cells, :); 
    end
end

%%
hold on
def_col = get(gca, 'ColorOrder');
for exp = experiment;
    %plot(Conv_Blocks, BPS{exp, 1}, 'Color', def_col(1,:))
    plot(Conv_Blocks, BPS{exp, 2}, 'Color', def_col(exp,:))
end

xlim([0 57])
% ylim([-1 1])
title('Fit Rank 2, No Coupling')
xlabel('Number of Blocks')
