% clear
datapath='/Volumes/Lab/Users/Nora/NSEM/GLM_Output/fixedSP_rk1_linear_MU_PS_CP_p8IDp8/standardparams_orig/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{2}='NSEM_mapPRJ/';
fittypepath{1}='WN_mapPRJ/';

for fittype=1
    for exp=1
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) '*.mat']);
        n_cells=length(matfiles);
        
        %Initialize matrices
%         mu=zeros(n_cells,1);
%         bps=zeros(n_cells,1);
%         type=zeros(n_cells,1);
%         PS=zeros(n_cells,120);
%         K_space=zeros(n_cells, 13, 13);
        K_time=zeros(n_cells, 30);
        
        % Collect info from files
        for file=1:n_cells
            disp(file)
            load([datapath fittypepath{fittype} exp_names(exp,:) matfiles(file).name]);
            % K_space(file,:,:)=fittedGLM.linearfilters.Stimulus.space_rk1;
            K_time(file,:)=fittedGLM.linearfilters.Stimulus.time_rk1;
%             mu(file)=fittedGLM.linearfilters.TonicDrive.Filter;
%             PS(file,:)=fittedGLM.linearfilters.PostSpike.Filter;
%             bps(file)=fittedGLM.xvalperformance.glm_normedbits;
%             if strcmp(matfiles(file).name(2), 'N')
%                 type(file)=1;
%             else
%                 type(file)=2;
%             end
        end
    end
end

%%

params = fittedGLM.rawfit.opt_params;
t1_idx = fittedGLM.rawfit.paramind.time1;
t2_idx = fittedGLM.rawfit.paramind.time2;
s1_idx = fittedGLM.rawfit.paramind.space1;
s2_idx = fittedGLM.rawfit.paramind.space2;
time = (1:length(t1_idx))*1/120*1000;
figure('Position', [100 100 700 200]);
subplot(1,3,1)
plot(time, params(t1_idx))
hold on
plot(time, params(t2_idx))
% title(cell_savename)
xlim([0 250])
xlabel('Time (ms)')
subplot(1,3,2)
imagesc(reshape(params(s1_idx),11,11))
caxis([0 0.5])
axis image
subplot(1,3,3)
imagesc(reshape(params(s2_idx),11,11))
caxis([0 0.5])
axis image
% print(['Users/Nora/Desktop/' cell_savename '.eps'], 'eps')

%%

plot(WN/max(WN), '--', 'LineWidth', 2)
hold on
plot(NSEM/max(NSEM), '--', 'LineWidth', 2)
plot(-params(t2_idx)/max(-params(t2_idx)), 'k')
plot(-params(t1_idx)/max(-params(t1_idx)), 'k')
legend('Mean Rank 1 WN Filter', 'Mean Rank 1 NSEM Filter', 'Sample Rank 2 NSEM Filters')
xlabel('Frame')
title('Multiple Timescales in the Temporal Part of the Stimulus Filter')
