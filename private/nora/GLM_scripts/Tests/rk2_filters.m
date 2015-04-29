% clear
datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_noCP_p8IDp8_shortlist/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{2}='NSEM_mapPRJ/';
fittypepath{1}='WN_mapPRJ/';

BPS_rk2 = zeros(100, 2);
plottin = 1;

for fittype=1
    for exp=1:4
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) '*.mat']);
        n_cells=length(matfiles);
        
        % Collect info from files
        for file=1:n_cells
            disp(file)
            load([datapath fittypepath{fittype} exp_names(exp,:) matfiles(file).name]);
            
            BPS_rk2(file, fittype) = fittedGLM.xvalperformance.glm_normedbits;
            
            % plotting
            if plottin
                time = (1:length(fittedGLM.linearfilters.Stimulus.time_rk1))*1/120*1000;
                figure('Position', [100 100 700 200]);
                subplot(1,3,1)
                plot(time, fittedGLM.linearfilters.Stimulus.time_rk1)
                hold on
                plot(time, fittedGLM.linearfilters.Stimulus.time_rk2)
                ylim([-1 1])
                title(fittedGLM.cellinfo.cell_savename)
                xlim([0 250])
                xlabel('Time (ms)')
                subplot(1,3,2)
                imagesc(fittedGLM.linearfilters.Stimulus.space_rk1)
                caxis([0 0.15])
                axis image
                subplot(1,3,3)
                imagesc(fittedGLM.linearfilters.Stimulus.space_rk2)
                caxis([0 0.15])
                axis image
                set(gcf,'PaperOrientation','landscape');
                set(gcf,'PaperPositionMode','auto')
                print('-dpsc',['/Users/Nora/Desktop/rk2/' exp_names(exp,:) fittedGLM.cellinfo.cell_savename '.ps'])
            end
            
        end
    end
end