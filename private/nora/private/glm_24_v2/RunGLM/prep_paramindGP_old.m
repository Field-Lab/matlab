% 2014-03-05 akheitman
% hack filled .. bt just make it a function to clear room in the main GLM
% code

function [paramind ] =  prep_paramindGP(GLMType, GLMPars)

numParams = 0;




if GLMType.CONVEX
    if GLMType.TonicDrive
        paramind.MU = 1;
        numParams = numParams+1;
    end
    if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')
        Xstart     = numParams + 1;  Xend = numParams +  GLMPars.stimfilter.frames;
        paramind.X = [Xstart:Xend];
        numParams  = numParams + GLMPars.stimfilter.frames;
    else
        error('you need to properly specifiy the stimfilter in prep_paramind')
    end

    if GLMType.PostSpikeFilter
        PSstart = numParams + 1;  PSend = numParams + GLMPars.spikefilters.ps.filternumber;
        paramind.PS = [PSstart  : PSend];
        numParams = numParams + GLMPars.spikefilters.ps.filternumber;
    end
    
elseif ~GLMType.CONVEX
    convParams = 0;
    if GLMType.TonicDrive
        paramind.convex.MU = 1;
        convParams = convParams+1;
    end
    if GLMType.PostSpikeFilter
        PSstart = convParams + 1;  PSend = convParams + GLMPars.spikefilters.ps.filternumber;
        paramind.convex.PS = [PSstart  : PSend];
        convParams = convParams + GLMPars.spikefilters.ps.filternumber;
    end
        
    if strcmp(GLMType.stimfilter_mode, 'rk1')
        Xstart = convParams + 1;  
        Xend   = convParams + (GLMPars.stimfilter.ROI_length^2) + GLMPars.stimfilter.frames;
        
        paramind.nonconvex.X      = [Xstart:Xend];
        paramind.nonconvex.space1 = [Xstart: ((Xstart-1) + (GLMPars.stimfilter.ROI_length^2))];
        paramind.nonconvex.time1  = [(Xstart + GLMPars.stimfilter.ROI_length^2) : Xend ];
        numParams = convParams  +  GLMPars.stimfilter.ROI_length^2 + GLMPars.stimfilter.frames;
    end
    
    paramind.convParams = convParams;
end

paramind.paramcount = numParams;


% Old way
%{ 
    if strcmp(Basepars.k_filtermode , 'raw')
        numParams  = 1 + Basepars.k_stimframes * (Basepars.k_spacepixels) + Basepars.ps_filternumber + length(Basepars.cp_Neighbors)*Basepars.cp_filternumber; 
        MU      = 1;   
        Lstart  = 2;               Lend  = (Lstart-1)  + Basepars.k_stimframes * (Basepars.k_spacepixels);
        PSstart = Lend  + 1;       PSend = (PSstart-1) + Basepars.ps_filternumber;
        CPstart = PSend + 1;       CPend = (CPstart-1) + length(Basepars.cp_Neighbors)*Basepars.cp_filternumber;  
        paramind.MU = MU;                MUind  = MU;
        paramind.L  = [Lstart  :  Lend]; Lind   = [Lstart  :  Lend];
        paramind.PS = [PSstart : PSend]; PSind  = [PSstart : PSend];
        paramind.CP = [CPstart : CPend]; CPind  = [CPstart : CPend];   
    elseif strcmp(Basepars.k_filtermode, 'rk2')
        numParams  = 1 + 2*( Basepars.k_stimframes + Basepars.k_spacepixels) + Basepars.ps_filternumber + length(Basepars.cp_Neighbors)*Basepars.cp_filternumber; 
        MU      = 1;   
        Lstart  = 2;                                    Lend  = (Lstart-1)  + 2*( Basepars.k_stimframes + Basepars.k_spacepixels);
        PSstart = Lend  + 1;                            PSend = (PSstart-1) + Basepars.ps_filternumber;
        CPstart = PSend + 1;                            CPend = (CPstart-1) + length(Basepars.cp_Neighbors)*Basepars.cp_filternumber;  
        paramind.MU = MU;                     MUind  = MU;
        paramind.L  = [Lstart  :  Lend];      Lind   = [Lstart  :  Lend];
        paramind.PS = [PSstart : PSend];      PSind  = [PSstart : PSend];
        paramind.CP = [CPstart : CPend];      CPind  = [CPstart : CPend]; 
        
        spixels = Basepars.k_spacepixels;     tframes = Basepars.k_stimframes;
        space1_start = Lstart;              space1_end = (space1_start-1) + spixels;
        time1_start  = space1_end + 1 ;      time1_end = (time1_start -1) + tframes;
        space2_start = time1_end+1;         space2_end = (space2_start-1) + spixels;
        time2_start  = space2_end + 1 ;      time2_end = (time2_start -1) + tframes;
        
        paramind.SPACE1 = [space1_start : space1_end];  SP1ind   = [space1_start : space1_end];
        paramind.TIME1  = [time1_start  :  time1_end];   T1ind   = [time1_start  :  time1_end];
        paramind.SPACE2 = [space2_start : space2_end];  SP2ind   = [space2_start : space2_end];
        paramind.TIME2  = [time2_start  :  time2_end];   T2ind   = [time2_start  :  time2_end];
	elseif strcmp(Basepars.k_filtermode , 'fixedSP' )
        numParams  = 1 + Basepars.k_stimframes + Basepars.ps_filternumber + length(Basepars.cp_Neighbors)*Basepars.cp_filternumber; 
        MU      = 1;   
        Lstart  = 2;                                    Lend  = (Lstart-1)  + Basepars.k_stimframes;
        PSstart = Lend  + 1;                            PSend = (PSstart-1) + Basepars.ps_filternumber;
        CPstart = PSend + 1;                            CPend = (CPstart-1) + length(Basepars.cp_Neighbors)*Basepars.cp_filternumber;  
        paramind.MU = MU;                     MUind  = MU;
        paramind.L  = [Lstart  :  Lend];      Lind   = [Lstart  :  Lend];
        paramind.PS = [PSstart : PSend];      PSind  = [PSstart : PSend];
        paramind.CP = [CPstart : CPend];      CPind  = [CPstart : CPend];         
        paramind.TIME1  = paramind.L;
        
	elseif strcmp(Basepars.k_filtermode , 'OnOff_hardrect_fixedSP_STA' )
        numParams  = 1 + 2* Basepars.k_stimframes + Basepars.ps_filternumber + length(Basepars.cp_Neighbors)*Basepars.cp_filternumber; 
        MU      = 1;   
        Lstart  = 2;                                    Lend  = (Lstart-1)  + 2* Basepars.k_stimframes;        
        PSstart = Lend  + 1;                            PSend = (PSstart-1) + Basepars.ps_filternumber;
        CPstart = PSend + 1;                            CPend = (CPstart-1) + length(Basepars.cp_Neighbors)*Basepars.cp_filternumber;  
        paramind.MU = MU;                     MUind  = MU;
        paramind.L  = [Lstart  :  Lend];      Lind   = [Lstart  :  Lend];
        paramind.PS = [PSstart : PSend];      PSind  = [PSstart : PSend];
        paramind.CP = [CPstart : CPend];      CPind  = [CPstart : CPend]; 
        
        paramind.TIME1 = [Lstart: ( (Lstart-1) + Basepars.k_stimframes) ];
        paramind.TIME2 = [((Lstart + Basepars.k_stimframes)) : Lend];
        
	elseif strcmp(Basepars.k_filtermode, 'STA')
        numParams  = 1 + 1 + Basepars.ps_filternumber + length(Basepars.cp_Neighbors)*Basepars.cp_filternumber; 
        MU      = 1;   
        Lscale  = 2;                     
        PSstart = Lscale + 1;                           PSend  = (PSstart-1) + Basepars.ps_filternumber;
        CPstart = PSend + 1;                            CPend  =  (CPstart-1) + length(Basepars.cp_Neighbors)*Basepars.cp_filternumber;  
        paramind.MU = MU;                     MUind  = MU;
        paramind.L  = Lscale;                 Lind   = Lscale;
        paramind.PS = [PSstart : PSend];      PSind  = [PSstart : PSend];
        paramind.CP = [CPstart : CPend];      CPind  = [CPstart : CPend]; 
    end
    if ~Basepars.Coupling
        numParams = numParams - length(Basepars.cp_Neighbors)*Basepars.cp_filternumber; 
        paramind = rmfield(paramind , 'CP' );
    end
    %{
    if Basepars.BiDirect_CP && Basepars.Coupling
         CPstart = CPstart;    CPend0 = CPend;
         cpparamsnum0 = 1+ CPend0 - CPstart;
         cpparamsnum  = 2 * cpparamsnum0;
         CPend        = CPstart + cpparamsnum - 1;
         paramind.CP = [CPstart : CPend ];
         numParams  =  cpparamsnum0 + numParams;
         CPind = [CPstart: CPend];
    end
    %}
    if Basepars.ps_FIX
        adjust                = Basepars.ps_filternumber - 1;
        paramind.PS = paramind.PS(1);
        PSind                 = paramind.PS(1);
        paramind.CP = paramind.CP - adjust;
        CPind                 = paramind.CP;
        numParams             = numParams - adjust;
    end    
    if Basepars.ps_FIX
        Basepars.ps_Filter  = zeros(Basepars.ps_filternumber,1);
    end
   %} 
   
    