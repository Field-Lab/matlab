% clear
datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/fixedSP_rk1_linear_MU_PS_CP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{2}='NSEM_mapPRJ/';
fittypepath{1}='WN_mapPRJ/';

count = 0;

for fittype=2
    for exp=1:4
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) '*.mat']);
        n_cells=length(matfiles);
        
        % CP=zeros(n_cells*6, 120);
        
        % Collect info from files
        for file=1:n_cells
            disp(file)
            load([datapath fittypepath{fittype} exp_names(exp,:) matfiles(file).name]);
            for pair = 1:6
                cp_max = max(fittedGLM.linearfilters.Coupling.Filter{pair});
                if  cp_max > 0.5
                    count = count + 1;
                    CP(count, :) = fittedGLM.linearfilters.Coupling.Filter{pair};
                end
            end
        end
    end
end

%%

plot(CP')
title('Normalized Coupling Filters')

%% PCA
[coeff, score, ~] = pca(CP);
waveform_NSEM = coeff(:,1);
init = score(:,1);
