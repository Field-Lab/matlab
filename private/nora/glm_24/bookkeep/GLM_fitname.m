% AKHeitman 2014-04-06
% concise naming scheme parameter folder name associated with a fit
% built in mechanism to locate older fits by optional condition
% Input the standard GLMType which dictate the fit_type


function GLM_fitname = GLM_fitname(GLMType , condition)
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
        core = sprintf('%s_SU'  , core);
    end
    if GLMType.CBP
        core = sprintf('%s_CBP'  , core);
    end
    if GLMType.Saccades
        core = sprintf('%s_S'  , core);
    end
    
    GLM_fitname_core = sprintf('%s_%s', core, GLMType.cone_sname);
    
    % modify core name is we are running a reduced debug mode
    if GLMType.debug  
        GLM_fitname_core = sprintf('debug_%s', GLM_fitname_core);
    end
    
    if isfield(GLMType, 'input_pt_nonlinearity') && GLMType.input_pt_nonlinearity
        GLM_fitname_core = sprintf('%s_stimnonlin_%s',GLM_fitname_core, GLMType.input_pt_nonlinearity_type)
    end
    
    % 
    if ~GLMType.specialchange
        GLM_fitname = sprintf('%s/standardparams',GLM_fitname_core)
    end
    if GLMType.specialchange && ~isempty(GLMType.specialchange_name)
        GLM_fitname = sprintf('%s/ChangeParams_%s',GLM_fitname_core,GLMType.specialchange_name )
        
    end
    
    
        
    
    
end

if exist('condition', 'var') && strcmp(condition, 'netappsnle')
    if  GLMType.Coupling,  Main_SolPars = sprintf('%s_ps%d_%s',GLMType.stimfilter_mode,GLMType.n_psf,'cpON');  end
    if ~GLMType.Coupling,  Main_SolPars = sprintf('%s_ps%d_%s',GLMType.stimfilter_mode,GLMType.n_psf,'cpOFF'); end

     if isfield(GLMType, 'rect') && GLMType.rect
        Main_SolPars = sprintf('%s_%s' , GLMType.rect_type , Main_SolPars );
    end 
    
	Other_SolPars = sprintf( 'bin%d_blk%d_tolfun%d',GLMType.binning,length(SPars.FitBlocks),GLMType.tolfun);
    if GLMType.Coupling, Other_SolPars = sprintf('%s_cpbasisno%d_neighborpertype%d', Other_SolPars,GLMType.n_cp, GLMType.neighbornumber); end
    if GLMType.debug, Other_SolPars = sprintf('debug_%s', Other_SolPars); end
    
    GLM_fitname = sprintf('%s/%s/%s', Main_SolPars, Other_SolPars,GLMType.cone_model);

end

end
