contrast_size = 11;
stimtype = 'NSEM';
exp = 2;
celltype = 'O';
datapath_NSEM='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_CP_p8IDp8/standardparams/NSEM_mapPRJ/';
datapath_WN='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_CP_p8IDp8/standardparams/WN_mapPRJ/';
exp_names=['2012-08-09-3/conv_blocks_57/';'2012-09-27-3/conv_blocks_57/';'2013-08-19-6/conv_blocks_57/'];

% Get all fits of that cell type
matfiles=dir([datapath_WN exp_names(exp,:) celltype '*.mat']);
n_cells=length(matfiles);

% Load one fittedGLM to get the experiment info
load([datapath exp_names(exp,:) matfiles(1).name]);
GLMType = fittedGLM.GLMType;
exp_nm = fittedGLM.cellinfo.exp_nm;
% Load and process stimulus
[StimulusPars, ~] = StimulusParams(exp_nm, stimtype, GLMType.map_type);
[movie, inputstats, ~]          = loadmoviematfile(exp_nm , stimtype, GLMType.cone_model,'testmovie');
testmovie             = movie{1}.matrix(:,:,StimulusPars.slv.testframes);

% Loop through the cells and find the residual spikes for each cell.
%error_WN = zeros(n_cells, 11400);
error_NSEM = zeros(n_cells, 34800);
%type = zeros(n_cells, 1);
%loc = zeros(n_cells, 2);
contrast_NSEM = zeros(n_cells, 34800);
for i_cell = 1:n_cells
    disp(i_cell)
    % Load fittedGLM
    load([datapath_NSEM exp_names(exp,:) matfiles(i_cell).name]);
    error_NSEM(i_cell, :) = log(glm_error(fittedGLM));
    ROI = ROI_coord(contrast_size, fittedGLM.cellinfo.slave_centercoord, StimulusPars.slv);
    contrast_temp = mean(reshape(double(testmovie(ROI.xvals, ROI.yvals, :)), [contrast_size^2, 3480]));
    contrast_NSEM(i_cell,:) = imresize(contrast_temp, [1, 34800],'nearest');
%     try
%         load([datapath_NSEM exp_names(exp,:) matfiles(i_cell).name]);
%         error_NSEM(i_cell, :) = glm_error(fittedGLM);
%         type(i_cell) = (fittedGLM.cell_savename(2) == 'N');
%         loc(i_cell,:) = [fittedGLM.cellinfo.slave_centercoord.x_coord, fittedGLM.cellinfo.slave_centercoord.y_coord];
%     catch
%    end
    
end


%%
% 
% if 0
% for i = 1:43
%     for j = 1:i
%         disp(i)
%         disp(j)
%         plot(error(i,:), error(j,:))
%         pause()
%     end
% end
% end

% %%
% centers = -200:200;
% counts_NSEM = hist(error_NSEM(2,:), centers);
% counts_WN = hist(error_WN(2,:), centers);
% plot(centers, counts_NSEM/34800)
% hold on
% plot(centers, counts_WN/11400)

% %% 
% % Try to find linear relationship
dur = 20000;
filter_L = 1000;
start_time = filter_L+100;
time = (start_time+1):(start_time+dur);
X = zeros(dur, filter_L);
for i = time
    X(i-start_time,:) = contrast_NSEM(2,(i-filter_L+1):i);
end
b = linsolve(X,error_NSEM(2,time)');
plot(b);


% % %% What about reverse corr?
% cvell = 1;
% filter_L = 1000;
% ETA = zeros(1,filter_L);
% control = zeros(1,filter_L);
% % contrast_NSEM = contrast_NSEM - mean(contrast_NSEM(1,:));
% for i = (filter_L+1):30000
%     ETA = error_NSEM(cvell,i) * contrast_NSEM(cvell,(i-filter_L+1):i) + ETA;
%     control =  contrast_NSEM(cvell,(i-filter_L+1):i) + control;
% end
% plot(ETA)
% figure;
% plot(control, '--')



