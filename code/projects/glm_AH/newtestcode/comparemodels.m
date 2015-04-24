% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS
% Modification of comparefit_template  started 2014-06-18


clear ; close all;
%
%Type = 'Stim_Nonlinearity'; modification = 'piece_linear_aboutmean'; doubleopt = true;
print_individual_comparison = false;
for i_fit = 1:2
    if i_fit == 1, GLMType.fit_type = 'WN';   GLMType.map_type = 'mapPRJ';end
    if i_fit == 2, GLMType.fit_type = 'NSEM'; GLMType.map_type = 'mapPRJ';end


%Type = 'Special_BiDirectional_StimFilter'; modification = 'linear_timekernel_shift0'; modification2 ='timekernel'; doubleopt = false;
%Type = 'BiDirectional_StimFilter'; modification = Type;  doubleopt = false;
%Type = 'ConeModelComparison';  modification2 ='timekernel' ; modification = 'linear_timekernel_shift0'; doubleopt = false;
%Type = 'TimeKernelComparison'; modification2 = 'timekernel'; modification = 'linear_timekernel_shift0'; doubleopt = false;
%Type = 'TimeKernelComparison'; modification2 = 'timekernelCONEMODEL'; modification = 'DimFlash_092413Fc12_shift0'; doubleopt = false;
%Type = 'Stim_Nonlinearity'; modification = 'log'; doubleopt = false;
%Type = 'Stim_Nonlinearity'; modification = 'exp'; doubleopt = false;
%Type = 'Stim_Nonlinearity'; modification = 'polynomial_androot_order2'; doubleopt = true;
%Type = 'Stim_Nonlinearity'; modification =
%'polynomial_androot_order2_search2'; doubleopt = true;
%Type = 'Stim_Nonlinearity'; modification = 'raisepower_meanafter'; doubleopt =true;
%Type = 'Stim_Nonlinearity'; modification = 'polynomial_order5_part4';
%doubleopt =true; doubleoptmanual = true; doubleopt_dim = 5; doubleoptcomp = false; %original = 'piece_linear_aboutmean';
%Type = 'Stim_Nonlinearity'; modification = 'piecelinear_fourpiece_eightlevels'; doubleopt =true; doubleoptmanual = true; doubleopt_dim = 5; doubleoptcomp = false
%Type = 'Stim_Nonlinearity'; modification =
%'oddfunc_powerraise_aboutmean';doubleopt =true; doubleopt_dim = 1;
%Type = 'Stim_Nonlinearity'; modification = 'piece_linear_shiftmean';original = 'piece_linear_aboutmean';doubleopt =true; doubleoptmanual = true; doubleopt_dim = 2;
%Type = 'Search_Method'; modification = 'inputNL'; nonlinearity = 'piece_linear_shiftmean'; oldsearch = 'standard'; newsearch = 'manual_grid'; doubleopt = false; doubleoptcomp = true;
%Type = 'FilterOutput_Nonlinearity'; modification ='piece_linear_aboutmean'; doubleopt = true;
%Type = 'FilterOutput_Nonlinearity'; modification ='oddfunc_powerraise_aboutmean';  doubleopt =true;
Type = 'GLM_StimFilter'; modification = 'rk2'; doubleopt = false; doubleoptcomp = false;
BD = NSEM_BaseDirectories;
baseoutput_dir = sprintf( '%s/Compare_Models', BD.NSEM_home);
exptests = [1 2 3 4];
i_exp = 1;  i_cell = 1;
cellselectiontype = 'shortlist';
%%%%%%%%%%%%%%%%%%%%%%%
% First set the GLMType
%%%%%%%%%%%%%%%%%%%%%%%
%GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
GLMType.cone_model = 'DimFlash_092413Fc12_shift0'; GLMType.cone_sname = 'timekernelCONEMODEL';
%GLMType.cone_model = 'linear_timekernel_shift0'; GLMType.cone_sname = 'timekernel';
%GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname = 'p8Mod1Max1e4p8';

%GLMType.fit_type = 'NSEM'; GLMType.map_type = 'mapPRJ';
%GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname =
%'p8Mod1Max1e4p8'; 
%GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
%GLMType.k_filtermode = 'OnOff_hardrect_fixedSP_STA'; GLMType.fixedSPlength = 13;  GLMType.fixedSP_nullpoint = 'mean'; 
GLMType.nullpoint = 'mean'; 
GLMType.debug = false; 
GLMType.specialchange = false;
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
%GLMType.stimfilter_mode = 'rk1';
GLMType.CONVEX = true;
GLMType.TonicDrive = true;
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.fixed_spatialfilter = true;
  

if exist('original', 'var') 
GLMType.input_pt_non_linearity = true;
GLMType.input_pt_non_linearity_type = original;
end
GLMType.fitname  = GLM_fitname(GLMType);


model_comparison.byexpnm        = cell(4,1);
model_comparison.comparisontype = Type; 
model_comparison.modification   = modification;
model_comparison.fit_type       = GLMType.fit_type;


%%%%%%%%%%%%%%%%%%%%%%%
% Set the comparison
%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Type, 'Search_Method')
    GLMType.DoubleOpt = true;
    if strcmp(modification, 'inputNL')
        GLMType.input_pt_nonlinearity      = true;
        GLMType.input_pt_nonlinearity_type = nonlinearity ;
        GLMType.fitname = GLM_fitname(GLMType);
        
        GLMType2 = GLMType;
        if strcmp(newsearch,'manual_grid')
            GLMType2.DoubleOpt_Manual = true;
            GLMType2.fitname = GLM_fitname(GLMType2);
        end
        
        model_comparison.xlabel         = 'regular fmincon search';
        model_comparison.ylabel         = 'grid search';
        model_comparison.title          = sprintf('Searching: %s vs %s', oldsearch, newsearch) 
        
        outputdir  = sprintf('%s/Stim_Nonlinearity/Search_Method/%s_v_%s', baseoutput_dir, oldsearch, newsearch);
        comparename = sprintf('InputSearch_Method');  
    end
end

if strcmp(Type, 'Stim_Nonlinearity')
    GLMType2     = GLMType;
    GLMType2.input_pt_nonlinearity      = true;
    GLMType2.input_pt_nonlinearity_type = modification ;
    
    if exist('doubleoptmanual','var') && doubleoptmanual
        GLMType2.DoubleOpt = true;
        GLMType2.DoubleOpt_Manual = true;
    end
    
    GLMType2.fitname                    = GLM_fitname(GLMType2); 
    
    model_comparison.xlabel         = 'regular stimulus';
    model_comparison.ylabel         = 'point nonlinearity';
    model_comparison.title          = sprintf('Point Nonlinearity: %s', modification) 
    
    outputdir  = sprintf('%s/Stim_Nonlinearity/%s/%s', baseoutput_dir, modification, GLMType.fitname);
    comparename = sprintf('%s: %s',  Type, modification);
end

if strcmp(Type, 'FilterOutput_Nonlinearity')
    GLMType2     = GLMType;
    GLMType2.postfilter_nonlinearity      = true;
    GLMType2.postfilter_nonlinearity_type = modification ;
    GLMType2.fitname                    = GLM_fitname(GLMType2); 
    outputdir  = sprintf('%s/FilterOutput_Nonlinearity/%s/%s', baseoutput_dir, modification, GLMType.fitname);
    
    model_comparison.xlabel         = 'standard model';
    model_comparison.ylabel         = 'output non-linearity';
    model_comparison.title          = sprintf('Filter Output Nonlinearity: %s', modification) 
    comparename = sprintf('%s: %s',  Type, modification);
end

if strcmp(Type, 'TimeKernelComparison')
    GLMType2     = GLMType;
    GLMType2.cone_model   = modification;
    GLMType2.cone_sname   = modification2; 
    GLMType2.fitname                    = GLM_fitname(GLMType2); 
    
    outputdir  = sprintf('%s/Stim_Filter/%s/%s', baseoutput_dir, modification, GLMType.fitname);
    model_comparison.xlabel         = sprintf('%s: %s', GLMType.fit_type, GLMType.cone_sname);
    model_comparison.ylabel         = sprintf('%s: %s', GLMType.fit_type, GLMType2.cone_sname);
    model_comparison.title          = sprintf('Comparing Cone Models') ;
    comparename = sprintf('%s: %s vs %s',  Type, GLMType.cone_sname,GLMType2.cone_sname);
    
end

if strcmp(Type, 'ConeModelComparison')
    GLMType2            = GLMType;
    GLMType2.cone_model = modification;
    GLMType2.cone_sname = modification2;
    GLMType2.fitname    = GLM_fitname(GLMType2);
    
    outputdir  = sprintf('%s/Stim_Filter/%s/%s', baseoutput_dir, modification, GLMType.fitname);
    model_comparison.xlabel         = sprintf('%s: %s', GLMType.fit_type, GLMType.cone_sname);
    model_comparison.ylabel         = sprintf('%s: %s', GLMType.fit_type, GLMType2.cone_sname);
    model_comparison.title          = sprintf('Comparing Cone Models') ;
    comparename = sprintf('%s: %s vs %s',  Type, GLMType.cone_sname,GLMType2.cone_sname);
end


if strcmp(Type,'BiDirectional_StimFilter')
    GLMType2                    = GLMType;
    GLMType2.specialchange      = true;
    GLMType2.specialchange_name = 'BiDirectional_StimFilter';
    GLMType2.fitname    = GLM_fitname(GLMType2);
    
    
    outputdir  = sprintf('%s/BiDirectional_StimFilter/%s/%s', baseoutput_dir, GLMType.fitname)
    model_comparison.xlabel         = sprintf('%s: %s', GLMType.fit_type, GLMType.cone_sname);
    model_comparison.ylabel         = sprintf('%s: %s with BiD Stim Filter', GLMType.fit_type, GLMType2.cone_sname);
    model_comparison.title          = sprintf('Bi-Directional Stim Filter') ;
    comparename = sprintf('%s',  Type);
end

if strcmp(Type, 'Special_BiDirectional_StimFilter')
    GLMType2                    = GLMType;
    GLMType2.specialchange      = true;
    GLMType2.specialchange_name = 'BiDirectional_StimFilter';
    GLMType2.cone_model = modification;
    GLMType2.cone_sname = modification2;
    GLMType2.fitname    = GLM_fitname(GLMType2);

    outputdir  = sprintf('%s/BiDirectional_StimFilter_andcmodel/%s/%s', baseoutput_dir, GLMType.fitname)
    model_comparison.xlabel         = sprintf('%s: %s', GLMType.fit_type, GLMType.cone_sname);
    model_comparison.ylabel         = sprintf('%s: %s with Bi-Directional Stim Filter', GLMType.fit_type, GLMType2.cone_sname);
    model_comparison.title          = sprintf('Bi-Directional Stim Filter and slow temporal filter') ;
    comparename = sprintf('%s',  Type);
end

if strcmp(Type, 'GLM_StimFilter')
    GLMType2                  = GLMType;
    GLMType2.stimfilter_mode  = modification;
    GLMType2.fitname          = GLM_fitname(GLMType2);
    outputdir                 = sprintf('%s/StimFilterChange/%s/%s', baseoutput_dir,modification,GLMType.fitname)
    model_comparison.xlabel         = sprintf('%s: %s', GLMType.fit_type, GLMType.stimfilter_mode);
    model_comparison.ylabel         = sprintf('%s: %s', GLMType.fit_type, GLMType2.stimfilter_mode);
    model_comparison.title          = sprintf('Varying GLM Stimulus Filter') ;
    comparename = sprintf('%s',  Type);
end


if ~exist(outputdir,'dir')
    if print_individual_comparison
        mkdir(outputdir);
    end
end

model_comparison.fitname1       = GLMType.fitname;
model_comparison.fitname2       = GLMType2.fitname;
model_comparison.fullGLMType    = GLMType;
model_comparison.note1          = 'rows are experiments, columns are the different models';





%% Load up parameters and xval scores across cells and experiements
for i_exp = exptests
    %% 
    expnumber = i_exp;
    [exp_nm,cells,expname]  = cell_list( expnumber, cellselectiontype);
    [StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
        
    rastdir = sprintf('%s/Rasters_%s/%s', outputdir,GLMType.fit_type,exp_nm) 
    if ~exist(rastdir,'dir'), mkdir(rastdir); end
    
    
    model_comparison.byexpnm{i_exp}.cells  = cells; 
    model_comparison.byexpnm{i_exp}.ONP    = zeros(1,length(cells));
    model_comparison.byexpnm{i_exp}.OFFP   = zeros(1,length(cells));
    model_comparison.byexpnm{i_exp}.normedbps_fit1             = zeros(1,length(cells));  
    model_comparison.byexpnm{i_exp}.normedbps_fit2             = zeros(1,length(cells));
    %%% if there is a parameter that we are optimizing %%%
    if exist('doubleopt','var') && doubleopt
        model_comparison.byexpnm{i_exp}.pt_nonlinearity_param      = cell(1,length(cells));
    end
    if exist('doubleoptcomp','var') && doubleoptcomp
        model_comparison.byexpnm{i_exp}.param_search1          = cell(1,length(cells));
        model_comparison.byexpnm{i_exp}.param_search2          = cell(1,length(cells));
    end
    model_comparison.byexpnm{i_exp}.Rast_Recorded              = cell(1,length(cells));
    model_comparison.byexpnm{i_exp}.Rast_fit1                  = cell(1,length(cells));
    model_comparison.byexpnm{i_exp}.Rast_fit2                  = cell(1,length(cells));
    
    %%%%%% Name and Create a Save Directory %%%%%%%%%%%
        
    inputs.exp_nm    = exp_nm; 
    inputs.map_type  = GLMType.map_type; 
    inputs.stim_type = GLMType.fit_type;
    inputs.fitname   = GLMType.fitname;
  
    d_save = NSEM_secondaryDirectories('savedir_GLMfit', inputs);  %clear inputs; 
    inputs2         = inputs;
    inputs2.fitname = GLMType2.fitname;
    d_save2         = NSEM_secondaryDirectories('savedir_GLMfit', inputs2);
    
   
 
    %% Load Cell Specific Elements   Spikes and STA
     for i_cell = 1:length(cells)
        cid = cells{i_cell};
        [celltype , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types);       
        
        eval(sprintf('load %s/%s.mat', d_save, cell_savename));
        fittedGLM1 = fittedGLM;
        
        eval(sprintf('load %s/%s.mat', d_save2, cell_savename));
        fittedGLM2 = fittedGLM;
        clear fittedGLM
        
        if strcmp(celltype, 'ON-Parasol')
            model_comparison.byexpnm{i_exp}.ONP(i_cell) = 1;
        elseif strcmp(celltype, 'OFF-Parasol')
            model_comparison.byexpnm{i_exp}.OFFP(i_cell) = 1;
        else
            error('messed up celltype naming');
        end
        model_comparison.byexpnm{i_exp}.normedbps_fit1(i_cell)          = fittedGLM1.xvalperformance.glm_normedbits;
        model_comparison.byexpnm{i_exp}.normedbps_fit2(i_cell)          = fittedGLM2.xvalperformance.glm_normedbits;
        
        %%% track the param if we are optimizing a param
        if doubleopt
            if strcmp(Type, 'Stim_Nonlinearity')
                model_comparison.byexpnm{i_exp}.pt_nonlinearity_param{i_cell}   = fittedGLM2.pt_nonlinearity_param;
            end
            if strcmp(Type, 'FilterOutput_Nonlinearity')
                if strcmp(modification,'piece_linear_aboutmean')
                    model_comparison.byexpnm{i_exp}.pt_nonlinearity_param{i_cell}   = fittedGLM2.GLMType.lcif_nonlinearity.increment_to_decrement;
                end
                if strcmp(modification,'oddfunc_powerraise_aboutmean')
                    model_comparison.byexpnm{i_exp}.pt_nonlinearity_param{i_cell}   = fittedGLM2.GLMType.lcif_nonlinearity.scalar_raisedpower;
                end
            end
        end
        
        if doubleoptcomp
            model_comparison.byexpnm{i_exp}.param_search1{i_cell}   = fittedGLM1.pt_nonlinearity_param;
            model_comparison.byexpnm{i_exp}.param_search2{i_cell}   = fittedGLM2.pt_nonlinearity_param;
        end
        if print_individual_comparison
            
            printname   = sprintf('%s/%s', rastdir,cell_savename);
            printmodelcomparison_withrate(fittedGLM1,fittedGLM2,comparename,printname);
        end

    end 
end
eval(sprintf('save %s/model_comparison_%s.mat model_comparison', outputdir, GLMType.fit_type));
%%  Plot change in xval score
    clf;
    subplot(3,1,1);
    axis off
    set(gca, 'fontsize', 12)
    c = 1;
    
    text(0, 1-0.1*c,sprintf('%s', model_comparison.title),'interpreter','none' )
    c = c+1;
    if ~strcmp(Type, 'ConeModelComparison')
    text(0, 1-0.1*c,sprintf('Fit by: %s, Cone Model: %s',GLMType.fit_type, GLMType.cone_model),'interpreter','none')
    else
        text(0, 1-0.1*c,sprintf('Fit by: %s',GLMType.fit_type, GLMType.cone_model),'interpreter','none')
    end
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit 1 on x-axis: %s',fittedGLM1.GLMType.fitname),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit 2 on y axis: %s',fittedGLM2.GLMType.fitname),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Dots are On-Parsols, Asterisks are Off-Parasols'))
    c = c+1;
    text(0, 1-0.1*c,sprintf('Each Color is a different experiment'))

    subplot(3,1,[2 3])
    plot(linspace(0,1,100), linspace(0,1,100),'k' ); hold on
    set(gca, 'fontsize', 14);
    set(gca,'xtick',[0:.2:.8]); set(gca,'ytick',[0:.2:.8]); 
    MS = 26;
    maxval = 0;
    minval = 1;
    for i_exp = 1:4
        fit1   = model_comparison.byexpnm{i_exp}.normedbps_fit1;
        fit2   = model_comparison.byexpnm{i_exp}.normedbps_fit2;  

        a = find(fit1<0); fit1(a) = 0; 
        b = find(  fit2<0);   fit2(b) = 0; 

        maxval= max( maxval , max(max(fit1),max(fit2)) );
        minval= min( minval , min(max(fit1),max(fit2)) );
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


    xlim([0,1.1*maxval]);
    ylim([0,1.1*maxval]); hold on;
    title('Normalized Bits per Spike')
    xlabel(model_comparison.xlabel, 'interpreter','none')
    ylabel(model_comparison.ylabel,'interpreter','none')    
    eval(sprintf('print -dpdf %s/XValScores_%s.pdf', outputdir,GLMType.fit_type));


%{
subplot(3,2,[4 6])
plot(linspace(0,1,100), linspace(0,1,100),'k' ); hold on
set(gca, 'fontsize', 14);
set(gca,'xtick',[0:.2:.8]); set(gca,'ytick',[0:.2:.8]); 
MS = 26;
maxval = 0;
minval = 0;
for i_exp = 1:4
    fit1   = model_comparison.byexpnm{i_exp}.normedbps_fit1;
    fit2   = model_comparison.byexpnm{i_exp}.normedbps_fit2;  
    
    a = find(fit1<=0);  
    b = find(  fit2<=0); 
    bad = union(a,b);
    
    param_vec = (fit2 - fit1) ./fit1 ;
    param_vec(bad) = 0;
    
    maxval = max(maxval, max(param_vec));
    minval = min(minval, min(param_vec));
    ONP  = find(model_comparison.byexpnm{i_exp}.ONP);
    OFFP = find(model_comparison.byexpnm{i_exp}.OFFP);
    
    
    
    
    ONP  = setdiff(ONP, bad);
    OFFP = setdiff(OFFP,bad);
    if i_exp == 1, plot(1*ones(1,length(ONP)), param_vec(ONP), 'r.','markersize',MS ); end
    if i_exp == 2, plot(2*ones(1,length(ONP)), param_vec(ONP), 'g.','markersize',MS ); end
    if i_exp == 3, plot(3*ones(1,length(ONP)), param_vec(ONP), 'b.','markersize',MS ); end
    if i_exp == 4, plot(4*ones(1,length(ONP)), param_vec(ONP), 'c.','markersize',MS ); end 
        
    if i_exp == 1, plot(1*ones(1,length(OFFP)), param_vec(OFFP), 'r*','markersize',MS ); end
    if i_exp == 2, plot(2*ones(1,length(OFFP)), param_vec(OFFP), 'g*','markersize',MS ); end
    if i_exp == 3, plot(3*ones(1,length(OFFP)), param_vec(OFFP), 'b*','markersize',MS ); end
    if i_exp == 4, plot(4*ones(1,length(OFFP)), param_vec(OFFP), 'c*','markersize',MS ); end     
end


xlim([0 5]);
ylim([minval,1.1*maxval]); hold on;
title('Percent Improvement in Normalized Bits Per Spike')
xlabel(model_comparison.xlabel, 'interpreter','none')
ylabel(model_comparison.ylabel,'interpreter','none')  
%}



%% show param clustering of 1-D double opt

if doubleopt  && (doubleopt_dim == 1)
    clf;

    subplot(3,1,1);
    axis off
    set(gca, 'fontsize', 12)
    c = 1;
    text(0, 1-0.1*c,sprintf('%s', model_comparison.title),'interpreter','none' )
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit by: %s, Cone Model: %s',GLMType.fit_type, GLMType.cone_model),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit name: %s',fittedGLM2.GLMType.fitname),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Dots are On-Parsols, Asterisks are Off-Parasols'))
    c = c+1;
    text(0, 1-0.1*c,sprintf('Each Color is a different experiment'))
    c = c+1;
    if strcmp(Type , 'Stim_Nonlinearity')
        text(0, 1-0.1*c,sprintf('Black line or is the parameter values at which we there is no stimulus modulation'))
    elseif strcmp(Type , 'FilterOutput_Nonlinearity') 
        text(0, 1-0.1*c,sprintf('Black line is the parameter value at which we there is no modulation of filter outputs'))
    end


    if strcmp(Type , 'Stim_Nonlinearity')
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'piece_linear_aboutmean')
            param_meaning = 'Ratio of rescaling of values above the mean to values below the mean';
            above_statement = 'Values ABOVE black line indicate comporatively more weight given to stim values ABOVE the original mean';
            below_statement = 'Values BELOW black line indicate comporatively more weight given to stim values BELOW the original mean';
        end
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'raisepower_meanafter')
            param_meaning = 'Power to which we raised the input stimulus (orig on 0 to 1 scale)';
            above_statement = 'Values ABOVE black line indicate comporatively more weight given to HIGHER stimulus values';
            below_statement = 'Values BELOW black line indicate comporatively more weight given to LOWER  stimulus values';
        end

        if strcmp(GLMType2.input_pt_nonlinearity_type, 'oddfunc_powerraise_aboutmean')
            param_meaning = 'Power to which we raised the input stimulus about mean made odd ';
            above_statement = 'Values ABOVE black line indicate we give comparatively more weight to extreme highs and extreme lows of the stimulus values';
            below_statement = 'Values BELOW black line indicate we give comparatively more weight to middle values of the stimulus value';
        end
    elseif strcmp(Type , 'FilterOutput_Nonlinearity')
        if strcmp(GLMType2.postfilter_nonlinearity_type, 'piece_linear_aboutmean')
            param_meaning = 'Ratio of rescaling of values of stimulus driven log CIF which are positive and negative';
            above_statement = 'Values ABOVE black line indicate comporatively more weight given to stim values ABOVE the original mean';
            below_statement = 'Values BELOW black line indicate comporatively more weight given to stim values BELOW the original mean';
        end
        if strcmp(GLMType2.postfilter_nonlinearity_type, 'oddfunc_powerraise_aboutmean')
            param_meaning = 'Ratio of rescaling of values of stimulus driven log CIF which are positive and negative';
            above_statement = 'Values ABOVE black line give added weight to extreme values (pos and neg) in lcif';
            below_statement = 'Values BELOW black line indicate gives added weight to values near 0';
        end    
    end

    c = c+1;
    text(0, 1-0.1*c,sprintf('Nonlinearity Parameter (y-axis): %s', param_meaning))



    subplot(3,1,[2 3])
    set(gca, 'fontsize', 14); hold on


    xlim0 = 0; xlim1 = 5;
    if strcmp(Type, 'Stim_Nonlinearity')
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'piece_linear_aboutmean')
            linear_point = 1;
        end
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'raisepower_meanafter')
            linear_point = 1;
        end
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'oddfunc_powerraise_aboutmean')
            linear_point = 1;
        end
    end

    if strcmp(Type, 'FilterOutput_Nonlinearity')
        if strcmp(GLMType2.postfilter_nonlinearity_type, 'piece_linear_aboutmean')
            linear_point = 1;
        end
        if strcmp(GLMType2.postfilter_nonlinearity_type, 'oddfunc_powerraise_aboutmean')
            linear_point = 1;
        end
    end


    LW = 2;    
    plot(linspace(xlim0,xlim1,100), linear_point*ones(1,100), 'k', 'linewidth', LW);
    xlim([xlim0,xlim1]);



    MS = 26;
    maxval = 0;
    minval = 1;
    for i_exp = 1:4
        param_vec = cell2mat([model_comparison.byexpnm{i_exp}.pt_nonlinearity_param(:)])
        maxval= max( maxval , max(param_vec) );
        minval= min( minval , min(param_vec) );

        ONP  = find(model_comparison.byexpnm{i_exp}.ONP);
        OFFP = find(model_comparison.byexpnm{i_exp}.OFFP);
        if i_exp == 1, plot(1*ones(1,length(ONP)), param_vec(ONP), 'r.','markersize',MS ); end
        if i_exp == 2, plot(2*ones(1,length(ONP)), param_vec(ONP), 'g.','markersize',MS ); end
        if i_exp == 3, plot(3*ones(1,length(ONP)), param_vec(ONP), 'b.','markersize',MS ); end
        if i_exp == 4, plot(4*ones(1,length(ONP)), param_vec(ONP), 'c.','markersize',MS ); end 

        if i_exp == 1, plot(1*ones(1,length(OFFP)), param_vec(OFFP), 'r*','markersize',MS ); end
        if i_exp == 2, plot(2*ones(1,length(OFFP)), param_vec(OFFP), 'g*','markersize',MS ); end
        if i_exp == 3, plot(3*ones(1,length(OFFP)), param_vec(OFFP), 'b*','markersize',MS ); end
        if i_exp == 4, plot(4*ones(1,length(OFFP)), param_vec(OFFP), 'c*','markersize',MS ); end          
    end
    set(gca,'xtick',[1 2 3 4]);
    set(gca,'xticklabel', {'expA', 'expB', 'expC', 'expD'})
    if strcmp(Type, 'Stim_Nonlinearity')
        title('Optimized Parameter of Stimulus Non-linearity')
    end
    if strcmp(Type, 'FilterOutput_Nonlinearity')
        title('Optimized Parameter of Filter Output Non-linearity')
    end
    ylabel('Non-linearity Parameter')    
    orient landscape
    eval(sprintf('print -dpdf %s/OptimizedParamClustering_%s.pdf', outputdir,GLMType.fit_type));
end


%% show param clustering for 2-D double opt
if doubleopt && doubleopt_dim == 2
    clf
    subplot(3,1,1);
    axis off
    set(gca, 'fontsize', 12)
    c = 1;
    text(0, 1-0.1*c,sprintf('%s', model_comparison.title),'interpreter','none' )
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit by: %s, Cone Model: %s',GLMType.fit_type, GLMType.cone_model),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit name: %s',fittedGLM2.GLMType.fitname),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Dots are On-Parsols, Asterisks are Off-Parasols'))
    c = c+1;
    text(0, 1-0.1*c,sprintf('Each Color is a different experiment'))
    c = c+1;
    if strcmp(Type , 'Stim_Nonlinearity')
        text(0, 1-0.1*c,sprintf('Black DOT is the parameter value with no stimulus modulation'))
    elseif strcmp(Type , 'FilterOutput_Nonlinearity') 
        text(0, 1-0.1*c,sprintf('Black DOT is the parameter value at which we there is no modulation of filter outputs'))
    end


    if strcmp(Type , 'Stim_Nonlinearity')
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'piece_linear_shiftmean')
            x_param_meaning = 'Ratio of rescaling of values above center to values below the mean';
            x_above_statement = 'Values ABOVE black line indicate comporatively more weight given to stim values ABOVE the original mean';
            x_below_statement = 'Values BELOW black line indicate comporatively more weight given to stim values BELOW the original mean';
            x_label           = 'Inc to Dec Ratio'
            
            y_param_meaning = 'Shift of Center away from additive mean';
            y_above_statement = 'Values ABOVE black line indicate comporatively more weight given to stim values ABOVE the original mean';
            y_below_statement = 'Values BELOW black line indicate comporatively more weight given to stim values BELOW the original mean';
            
            y_label          = 'New Center';
            
            no_modulation = [1,0];
            xlim0 = 0; xlim1 = 2;
            ylim0 = -.25; ylim1 = .25; 
        end  
    end

    c = c+1;
    text(0, 1-0.1*c,sprintf('Nonlinearity Parameter(x-axis): %s', x_param_meaning));
    c = c+1;
    text(0, 1-0.1*c,sprintf('Nonlinearity Parameter(y-axis): %s', y_param_meaning));
    
    subplot(3,1,[2 3]); xlim([xlim0,xlim1]); ylim([ylim0,ylim1]); hold on
    plot(no_modulation(1), no_modulation(2), 'k.', 'markersize',40); 
    plot(linspace(xlim0,xlim1,100), no_modulation(2)*ones(1,100),'k');
    plot(no_modulation(1)*ones(1,100),linspace(ylim0,ylim1,100), 'k');

    MS = 26;
    maxval = 0;
    minval = 1;
    for i_exp = 1:4
        param_vec = cell2mat([model_comparison.byexpnm{i_exp}.pt_nonlinearity_param(:)]);
        maxval= max( maxval , max(param_vec) );
        minval= min( minval , min(param_vec) );

        ONP  = find(model_comparison.byexpnm{i_exp}.ONP);
        OFFP = find(model_comparison.byexpnm{i_exp}.OFFP);
        if i_exp == 1, plot( param_vec(ONP,1), param_vec(ONP,2), 'r.','markersize',MS ); end
        if i_exp == 2, plot( param_vec(ONP,1), param_vec(ONP,2), 'g.','markersize',MS ); end
        if i_exp == 3, plot( param_vec(ONP,1), param_vec(ONP,2), 'b.','markersize',MS ); end
        if i_exp == 4, plot( param_vec(ONP,1), param_vec(ONP,2), 'c.','markersize',MS ); end 

        if i_exp == 1, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'r*','markersize',MS ); end
        if i_exp == 2, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'g*','markersize',MS ); end
        if i_exp == 3, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'b*','markersize',MS ); end
        if i_exp == 4, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'c*','markersize',MS ); end          
    end
    if strcmp(Type, 'Stim_Nonlinearity')
        title('Optimized Parameter of Stimulus Non-linearity')
    end
    if strcmp(Type, 'FilterOutput_Nonlinearity')
        title('Optimized Parameter of Filter Output Non-linearity')
    end
    ylabel(y_label); xlabel(x_label)    
    orient landscape
    eval(sprintf('print -dpdf %s/OptimizedParamClustering_%s.pdf', outputdir,GLMType.fit_type));
end
%%
if doubleoptcomp 
    clf
    subplot(3,1,1);
    axis off
    set(gca, 'fontsize', 12)
    c = 1;
    text(0, 1-0.1*c,sprintf('%s', model_comparison.title),'interpreter','none' )
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit by: %s, Cone Model: %s',GLMType.fit_type, GLMType.cone_model),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit name: %s',fittedGLM2.GLMType.fitname),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Dots are On-Parsols, Asterisks are Off-Parasols'))
    c = c+1;
    text(0, 1-0.1*c,sprintf('Each Color is a different experiment'))
    c = c+1;
    if strcmp(modification,'inputNL')
        text(0, 1-0.1*c,sprintf('Black DOT is the parameter value with no stimulus modulation'))
    elseif strcmp(Type , 'FilterOutput_Nonlinearity') 
        text(0, 1-0.1*c,sprintf('Black DOT is the parameter value at which we there is no modulation of filter outputs'))
    end


    if strcmp(modification,'inputNL')
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'piece_linear_shiftmean')
            x_param_meaning = 'Ratio of rescaling of values above center to values below the mean';
            x_above_statement = 'Values ABOVE black line indicate comporatively more weight given to stim values ABOVE the original mean';
            x_below_statement = 'Values BELOW black line indicate comporatively more weight given to stim values BELOW the original mean';
            x_label           = 'Inc to Dec Ratio'
            
            y_param_meaning = 'Shift of Center away from additive mean';
            y_above_statement = 'Values ABOVE black line indicate comporatively more weight given to stim values ABOVE the original mean';
            y_below_statement = 'Values BELOW black line indicate comporatively more weight given to stim values BELOW the original mean';
            
            y_label          = 'New Center';
            
            no_modulation = [0,0];
            xlim0 = -1; xlim1 = 1;
            ylim0 = -.3; ylim1 = .3; 
        end  
    end

    c = c+1;
    text(0, 1-0.1*c,sprintf('Nonlinearity Parameter(x-axis): %s', x_param_meaning));
    c = c+1;
    text(0, 1-0.1*c,sprintf('Nonlinearity Parameter(y-axis): %s', y_param_meaning));
    
    subplot(3,1,[2 3]); xlim([xlim0,xlim1]); ylim([ylim0,ylim1]); hold on
    %plot(no_modulation(1), no_modulation(2), 'k.', 'markersize',40); 
    plot(linspace(xlim0,xlim1,100), no_modulation(2)*ones(1,100),'k');
    plot(no_modulation(1)*ones(1,100),linspace(ylim0,ylim1,100), 'k');

    MS = 26;
    maxval = 0;
    minval = 1;
    for i_exp = 1:4
        param_vec = cell2mat([model_comparison.byexpnm{i_exp}.param_search2(:)]) - cell2mat([model_comparison.byexpnm{i_exp}.param_search1(:)])
      %  maxval= max( maxval , max(param_vec) );
       % minval= min( minval , min(param_vec) );

        ONP  = find(model_comparison.byexpnm{i_exp}.ONP);
        OFFP = find(model_comparison.byexpnm{i_exp}.OFFP);
        if i_exp == 1, plot( param_vec(ONP,1), param_vec(ONP,2), 'r.','markersize',MS ); end
        if i_exp == 2, plot( param_vec(ONP,1), param_vec(ONP,2), 'g.','markersize',MS ); end
        if i_exp == 3, plot( param_vec(ONP,1), param_vec(ONP,2), 'b.','markersize',MS ); end
        if i_exp == 4, plot( param_vec(ONP,1), param_vec(ONP,2), 'c.','markersize',MS ); end 

        if i_exp == 1, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'r*','markersize',MS ); end
        if i_exp == 2, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'g*','markersize',MS ); end
        if i_exp == 3, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'b*','markersize',MS ); end
        if i_exp == 4, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'c*','markersize',MS ); end          
    end
    title('Manual Params - FMINCON Params')
    ylabel(y_label); xlabel(x_label)    
    orient landscape
    eval(sprintf('print -dpdf %s/SearchDifferences_%s.pdf', outputdir,GLMType.fit_type));
end
end






