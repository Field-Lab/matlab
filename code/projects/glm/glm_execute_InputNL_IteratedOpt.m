function [fittedGLM] = glm_execute_InputNL_IteratedOpt(GLMType,fitspikes,fitmovie,testspikes_raster,testmovie,inputstats,glm_cellinfo,neighborspikes,troubleshoot)

GLMPars = GLMParams;
if strcmp(GLMType.input_pt_nonlinearity_type, 'powerraise')
    GLMPars.others.point_nonlinearity.powerraise = 1;
end

fittedGLM = [];



end
