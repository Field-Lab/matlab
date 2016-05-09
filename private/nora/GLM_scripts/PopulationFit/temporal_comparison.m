
clear

%%
datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/fixedSP_rk1_linear_MU_PS_CP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fit_type{1}='WN';
fit_type{2}='NSEM';
fittypepath{2}=[fit_type{2} '_mapPRJ/'];
fittypepath{1}=[fit_type{1} '_mapPRJ/'];

for Iexp = 2
    for fittype = 1
        count = 0;
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(Iexp,:) 'ON*.mat']);
        n_cells=length(matfiles);
        
        %Initialize matrices
        K_time=zeros(n_cells, 30);
        cid = zeros(n_cells,1);
        % Collect info from files
        for file=1:n_cells
            disp(file)
            load([datapath fittypepath{fittype} exp_names(Iexp,:) matfiles(file).name]);
            K_time(file,:)=fittedGLM.linearfilters.Stimulus.time_rk1;
            cid(file) = fittedGLM.cellinfo.cid;
        end
    end
end


%%
datarun = load_data('2012-09-27-3/data003');
datarun = load_params(datarun);
STA_time=zeros(n_cells, 30, 3);
for i =1:n_cells
    STA_time(i,:, 1) = datarun.vision.timecourses(datarun.cell_ids == cid(i)).r;
    STA_time(i,:, 2) = datarun.vision.timecourses(datarun.cell_ids == cid(i)).g;
    STA_time(i,:, 3) = datarun.vision.timecourses(datarun.cell_ids == cid(i)).b;
end
