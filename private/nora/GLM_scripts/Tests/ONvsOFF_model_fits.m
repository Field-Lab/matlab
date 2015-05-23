% clear
datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_noCP_p8IDp8_shortlist/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{2}='NSEM_mapPRJ/';
fittypepath{1}='WN_mapPRJ/';

for fittype=2
    
    for exp=1:4
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) '*.mat']);
        n_cells=length(matfiles);
        
        % Collect info from files
        for file=1:n_cells
            disp(file)
            load([datapath fittypepath{fittype} exp_names(exp,:) matfiles(file).name]);
            if strcmp(matfiles(file).name(2), 'N')
                type{exp}(file)=1;
            else
                type{exp}(file)=2;
            end
            BPS{exp, fittype}(file) = fittedGLM.xvalperformance.glm_normedbits;
        end
    end
end

%%

[countsON, centersON] = hist(BPS{2,2}(type{2} == 1));
[countsOFF, centersOFF] = hist(BPS{2,2}(type{2} == 2), centersON);

plot(centersON, countsON)
hold on
plot(centersOFF, countsOFF)
legend('ON', 'OFF')