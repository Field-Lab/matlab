% AKHeitman 2014-04-06
% concise naming scheme parameter folder name associated with a fit
% built in mechanism to locate older fits by optional condition
% Input the standard GLMType which dictate the fit_type


function GLM_fitname = GLM_fitname(GLMType)
% condition = 'netappsnle'  just a hack to load up old fits

if ~exist('condition','var') 
    
    % find the major model components define core fit name
    core = GLMType.stimfilter_mode;
    if GLMType.TonicDrive
        core = sprintf('%s_MU'  , core);
    else
        core = sprintf('%s_noMU', core);
    end    
    if GLMType.PostSpikeFilter
        core = sprintf('%s_PS'  , core);
    else
        core = sprintf('%s_noPS', core);
    end 
    if GLMType.CouplingFilters
        core = sprintf('%s_CP'  , core);
    else
        core = sprintf('%s_noCP', core);
    end
    if GLMType.Subunits
        core = sprintf('%s_SU%s'  , core, GLMType.Subunit_NL);
    end
    if GLMType.Contrast
        core = sprintf('%s_C'  , core); 
    end
    if isfield(GLMType, 'Saccades')
        core = sprintf('%s_SA'  , core); 
    end
    if GLMType.STA_init 
        core = sprintf('%s_init', core);
    end
    GLM_fitname_core = sprintf('%s_%s', core, GLMType.cone_sname);
    
    % modify core name is we are running a reduced debug mode
    if isfield(GLMType,'debug') && GLMType.debug  
        GLM_fitname_core = sprintf('debug_%s', GLM_fitname_core);
    end
    
    if isfield(GLMType, 'input_pt_nonlinearity') && GLMType.input_pt_nonlinearity
        GLM_fitname_core = sprintf('%s_stimnonlin_%s',GLM_fitname_core, GLMType.input_pt_nonlinearity_type);
    end
    
    if isfield(GLMType, 'postfilter_nonlinearity') && GLMType.postfilter_nonlinearity
        GLM_fitname_core = sprintf('%s_postfilterNL_%s',GLM_fitname_core, GLMType.postfilter_nonlinearity_type);
    end
    
    % NB time filter
    if isfield(GLMType,'timefilter')
        GLM_fitname_core = [GLM_fitname_core, GLMType.timefilter];
    end
    
    % See if there are special changes to the core parameters
    if ~isfield(GLMType, 'specialchange') || ~GLMType.specialchange
        GLM_fitname = sprintf('%s/standardparams',GLM_fitname_core);
    else
        GLM_fitname = sprintf('%s/ChangeParams_%s',GLM_fitname_core,GLMType.specialchange_name ) ;
    end
    if isfield(GLMType,'DoubleOpt_Manual') && GLMType.DoubleOpt_Manual
        GLM_fitname = sprintf('%s/Man_DoubleOpt_standardparams',GLM_fitname_core);
    end
    
    %
    if isfield(GLMType, 'name_change')
        GLM_fitname = [GLM_fitname GLMType.name_change];
    end
    

   
end


end