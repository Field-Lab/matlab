% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS
baseoutput_dir = '/Users/akheitman/NSEM_Home/CrossStim_Performance';

% SETUP cells and experiments, the TYPE of GLM (GLMType) 
BD = NSEM_BaseDirectories;
exptests = [1 2 3 4];
cellselectiontype = 'shortlist';%cellselectiontype = 'debug';
GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
%GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname = 'p8Mod1Max1e4p8';
%GLMType.k_filtermode = 'OnOff_hardrect_fixedSP_STA'; GLMType.fixedSPlength = 13;  GLMType.fixedSP_nullpoint = 'mean'; 
GLMType.nullpoint = 'mean'; 
GLMType.fit_type = 'NSEM'; GLMType.map_type = 'mapPRJ';
GLMType.debug = false; 
GLMType.specialchange = false;
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.CONVEX = false;
GLMType.TonicDrive = true;
GLMTYpe.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.fixed_spatialfilter = true;
GLMType.func_sname = 'glmwrap_23';
GLMType.fullmfilename =mfilename('fullpath'); 
i_exp = 1; i_cell = 1;

GLMType.fitname  = GLM_fitname(GLMType);  


GLMType2 = GLMType;
GLMType2.input_pt_nonlinearity      = true;
GLMType2.input_pt_nonlinearity_type = 'piece_linear_aboutmean';
GLMType2.fitname  = GLM_fitname(GLMType2);  

comparedir  = '/Users/akheitman/NSEM_Home/Compare_Models';
comparename = 'DO_piece_linear_aboutmean_NSEM';




model_comparison.byexpnm     = cell(4,1)
model_comparison.fit_type    = 'NSEM';
model_comparison.fitname1    = GLMType.fitname;
model_comparison.fitname2    = GLMType2.fitname;
model_comparison.xlabel      = 'no stim modulation';
model_comparison.ylabel      = 'point nonlinearity to stim'
model_comparison.title       = 'piecewise linear about mean';
model_comparison.fullGLMType = GLMType;
model_comparison.note1 = 'rows are experiments, columns are the different models'


savedir = sprintf('%s/%s', comparedir, comparename);
if ~exist(savedir); mkdir(savedir); end



%%
for i_exp = exptests
    %% 
    expnumber = i_exp;
    [exp_nm,cells,expname]  = cell_list( expnumber, cellselectiontype);
    if i_exp == 4, cells = cells(1:6); end
    [StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
    SPars = StimulusPars.slv;
    
    model_comparison.byexpnm{i_exp}.cells  = cells; 
    model_comparison.byexpnm{i_exp}.ONP    = zeros(1,length(cells));
    model_comparison.byexpnm{i_exp}.OFFP   = zeros(1,length(cells));
    model_comparison.byexpnm{i_exp}.normedbps_fit1             = zeros(1,length(cells));  
    model_comparison.byexpnm{i_exp}.normedbps_fit2             = zeros(1,length(cells));
    model_comparison.byexpnm{i_exp}.pt_nonlinearity_param      = cell(1,length(cells));      
    %%%%%% Name and Create a Save Directory %%%%%%%%%%%
        
    inputs.exp_nm    = exp_nm; 
    inputs.map_type  = GLMType.map_type; 
    inputs.stim_type = GLMType.fit_type;
    inputs.fitname   = GLMType.fitname;
  
    d_save = NSEM_secondaryDirectories('savedir_GLMfit', inputs);  %clear inputs; 
    
    inputs2         = inputs;
    inputs2.fitname = GLMType2.fitname;
    d_save2 = NSEM_secondaryDirectories('savedir_GLMfit', inputs2);
    
   
 
    %% Load Cell Specific Elements   Spikes and STA

    

    for i_cell = 1:length(cells)
        
        
        cid = cells{i_cell};
        [celltype , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types);       
        
        eval(sprintf('load %s/%s.mat', d_save, cell_savename));
        fittedGLM1 = fittedGLM;
        
        eval(sprintf('load %s/%s.mat', d_save2, cell_savename));
        fittedGLM2 = fittedGLM;

        if strcmp(celltype, 'ON-Parasol')
            model_comparison.byexpnm{i_exp}.ONP(i_cell) = 1;
        elseif strcmp(celltype, 'OFF-Parasol')
            model_comparison.byexpnm{i_exp}.OFFP(i_cell) = 1;
        else
            error('messed up celltype naming');
        end
        
        model_comparison.byexpnm{i_exp}.normedbps_fit1(i_cell)          = fittedGLM1.xvalperformance.glm_normedbits;
        model_comparison.byexpnm{i_exp}.normedbps_fit2(i_cell)          = fittedGLM2.xvalperformance.glm_normedbits;
        model_comparison.byexpnm{i_exp}.pt_nonlinearity_param{i_cell}   = fittedGLM2.pt_nonlinearity_param;
        
        
        
        %fittedGLM.xvalperformance.rasters  .recorded  .glm_sim   .bintime

    end
    
end
eval(sprintf('save %s/model_comparison.mat model_comparison', savedir));
%%
clf;

xlim([0,.8]);
ylim([0,.8]); hold on;
plot(linspace(0,1,100), linspace(0,1,100),'k' );
set(gca, 'fontsize', 14);
set(gca,'xtick',[0:.2:.8]); set(gca,'ytick',[0:.2:.8]); 
MS = 26;
for i_exp = 1:4
    fit1   = model_comparison.byexpnm{i_exp}.normedbps_fit1;
    fit2   = model_comparison.byexpnm{i_exp}.normedbps_fit2;  
    a = find(fit1<0); fit1(a) = 0; 
    b = find(  fit2<0);   fit2(b) = 0; 
    ONP  = find(model_comparison.byexpnm{i_exp}.ONP);
    OFFP = find(model_comparison.byexpnm{i_exp}.OFFP);
    if i_exp == 1, plot(fit1(ONP), fit2(ONP),'r.','markersize',MS ); end
    if i_exp == 2, plot(fit1(ONP), fit2(ONP),'g.','markersize',MS ); end
    if i_exp == 3, plot(fit1(ONP), fit2(ONP),'b.','markersize',MS ); end
    if i_exp == 4, plot(fit1(ONP), fit2(ONP),'c.','markersize',MS ); end 
    
    
    if i_exp == 1, plot(fit1(OFFP), fit2(OFFP),'r*','markersize',MS ); end
    if i_exp == 2, plot(fit1(OFFP), fit2(OFFP),'g*','markersize',MS ); end
    if i_exp == 3, plot(fit1(OFFP), fit2(OFFP),'b*','markersize',MS ); end
    if i_exp == 4, plot(fit1(OFFP), fit2(OFFP),'c*','markersize',MS ); end 
    
    
end
title(model_comparison.title)
xlabel(model_comparison.xlabel)
ylabel(model_comparison.ylabel)
    
orient landscape
eval(sprintf('print -dpdf %s/%s-Summary.pdf', savedir,comparename));



%eval(sprintf('print -dpdf WN_outperforms_NSEM.pdf'));
%}







