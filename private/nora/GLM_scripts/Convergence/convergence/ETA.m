datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{2}='NSEM_mapPRJ/';
fittypepath{1}='WN_mapPRJ/';

Conv_Blocks = [1 2 3 5 7 9 11 20 30 40 50 57];
conv_fit = 20; 

for fittype=1
    for exp=1
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) 'conv_blocks_20/*.mat']);
        
        % Collect info from files
        for file=1:2
            load([datapath fittypepath{fittype} exp_names(exp,:) 'conv_blocks_20/' matfiles(file).name]);
            
            spikes_home = '/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/';
            load([spikes_home exp_names(exp,:) fittypepath{fittype} 'organizedspikes_' matfiles(file).name], 'organizedspikes');
            [blockedmoviecell, inputstats, ~] = loadmoviematfile(exp_names(exp,1:(end-1)) , fittypepath{fittype}(1:(end-8)) ,'8pix_Identity_8pix', 'fitmovie');
            
            eta_blocks = (conv_fit+1):length(blockedmoviecell);
            spikes = organizedspikes.block.t_sp_withinblock(2:2:end);
            
            
            % I think we might need to record drive instead of comparing
            % two spike trains? 
            for i = 1%eta_blocks
                X=eval_xvalperformance_conv(fittedGLM, spikes{i}, blockedmoviecell{i}.matrix, inputstats, 0);
            end
        end
    end
end


%%
hold on
spk_tm = find(X.rasters.recorded);
plot(spk_tm, ones(length(spk_tm),1), '.')
spk_tm = find(X.rasters.glm_sim);
plot(spk_tm, 0.9*ones(length(spk_tm),1), '.')
ylim([0.5 1.5])
xlim([0 3000])

%%
con_length = 2000;
PSTH_rec = conv(X.rasters.recorded, ones(con_length,1));
PSTH_sim = conv(X.rasters.glm_sim, ones(con_length,1));
plot(PSTH_rec)
hold on
plot(PSTH_sim)