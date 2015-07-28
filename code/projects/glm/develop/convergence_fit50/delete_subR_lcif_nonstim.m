function [lcif_nonstim] = subR_lcif_nonstim(pstar, GLMType,GLMPars,fitspikes,t_bin,bins)
% AKHEITMAN 2015-07-15  extract non-stim components of lcif


% Construct PS_Basis
if GLMType.PostSpikeFilter
    basis_params  = GLMPars.spikefilters.ps;
    ps_basis      = prep_spikefilterbasisGP(basis_params,t_bin);
    

    if isfield(GLMType, 'special_arg') && isfield(GLMType.special_arg,'PS_Constrain')
            ps_basis_0 = ps_basis; clear ps_basis
            v        = sum(ps_basis_0,1);
            v        = v / norm(v) ;
            orthog_v = null(v);
            COB      = [v', orthog_v] ;
            ps_basis = (inv(COB) * ps_basis_0')' ;
    end

            
end

% Find PS_bin and MU_bin (pre-param multiply)
t_bin        = t_bin;
home_sptimes = fitspikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins  = home_spbins(find(home_spbins < bins) );
if GLMType.PostSpikeFilter
    basis         = ps_basis';
    PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
end
if GLMType.TonicDrive
    MU_bin = ones(1,bins);
end


% Multiply by param for final log-cif values
[paramind] =  prep_paramindGP(GLMType, GLMPars); 
pstar = pstar';
total_lcif = 0;
if isfield(paramind, 'MU')
    lcif_nonstim.components.mu   = pstar(paramind.MU)*MU_bin;
    total_lcif = total_lcif + lcif_nonstim.components.mu;
end
if isfield(paramind, 'PS')
    lcif_nonstim.components.ps   = pstar(paramind.PS)*PS_bin;
    total_lcif = total_lcif +  lcif_nonstim.components.ps;
end
lcif_nonstim.total = total_lcif;



end