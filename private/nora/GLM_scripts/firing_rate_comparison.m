
clear

%%
datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/fixedSP_rk1_linear_MU_PS_CP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fit_type{1}='WN';
fit_type{2}='NSEM';
fittypepath{2}=[fit_type{2} '_mapPRJ/'];
fittypepath{1}=[fit_type{1} '_mapPRJ/'];

for Iexp = 2
    for fittype = 2        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(Iexp,:) 'ON*.mat']);
        n_cells=length(matfiles);

        % Collect info from files
        for file=1:n_cells
            disp(file)
            load([datapath fittypepath{fittype} exp_names(Iexp,:) matfiles(file).name]);
            mu_on(file)=fittedGLM.linearfilters.TonicDrive.Filter;
        end
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(Iexp,:) 'OFF*.mat']);
        n_cells=length(matfiles);

        % Collect info from files
        for file=1:n_cells
            disp(file)
            load([datapath fittypepath{fittype} exp_names(Iexp,:) matfiles(file).name]);
            mu_off(file)=fittedGLM.linearfilters.TonicDrive.Filter;
        end
    end
end
