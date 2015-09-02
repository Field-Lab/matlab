clear
datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/fixedSP_rk1_linear_MU_PS_noCP_init_p8IDp8fit/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/'; '2013-10-10-0/'];
fittypepath{2}='NSEM_mapPRJ/';
fittypepath{1}='WN_mapPRJ/';
on_count = 0;
off_count = 0;

for fittype=2
    
    for exp=1
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) '*.mat']);
        n_cells=length(matfiles);
        
        % Collect info from files
        for file=1:n_cells
            if ~mod(file, 10); disp(file); end
            load([datapath fittypepath{fittype} exp_names(exp,:) matfiles(file).name]);
            if strcmp(matfiles(file).name(2), 'N')
                type{exp}(file)=1;
                on_count = on_count+1;
            else
                type{exp}(file)=2;
                off_count = off_count+1;
            end
            BPS{exp, fittype}(file) = fittedGLM.xvalperformance.glm_normedbits;
        end
    end
end

%%

%exp = 1; 
[countsOFF, centersOFF] = hist(BPS{exp,fittype}(type{exp} == 2), 15);
[countsON, centersON] = hist(BPS{exp,fittype}(type{exp} == 1), centersOFF);

plot(centersON, countsON/on_count, 'LineWidth', 2)
hold on
plot(centersOFF, countsOFF/off_count, 'LineWidth', 2)
legend('ON', 'OFF')
title('BPS Distribution')