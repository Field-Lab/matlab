
if ~GLMType.specialchange
    GLMPars = GLMParams;
end
if GLMType.specialchange
    GLMPars = GLMParams(GLMType.specialchange_name);
end

if GLMType.debug, GLMPars.tolfun = 2; end
if ~GLMType.Coupling
        GLMPars.n_cp = 0;
        GLMPars.BiDirect_Cp = false;
end