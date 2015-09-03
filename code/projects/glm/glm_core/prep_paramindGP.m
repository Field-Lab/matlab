% 2014-03-05 akheitman
% hack filled .. bt just make it a function to clear room in the main GLM
% code

function [paramind ] =  prep_paramindGP(GLMType, GLMPars)

numParams = 0;


% Assign Parameter Indices.. First to elements that won't brek convexity
% IE. Tonic Drive, Spike filters.. Whatever
if GLMType.TonicDrive
	paramind.MU = 1;
	numParams = numParams+1;
end
if GLMType.PostSpikeFilter
	PSstart = numParams + 1;  PSend = numParams + GLMPars.spikefilters.ps.filternumber;
	paramind.PS = [PSstart  : PSend];
	numParams = numParams + GLMPars.spikefilters.ps.filternumber;
end
% NBCoupling 2015-04-20
if GLMType.CouplingFilters
    for j_params=1:GLMPars.spikefilters.cp.n_couplings
        CPstart = numParams + 1;  CPend = numParams + GLMPars.spikefilters.cp.filternumber;
        paramind.CP{j_params} = [CPstart  : CPend]; % paramind.CP has numbers for each coupled cell
        numParams = numParams + GLMPars.spikefilters.cp.filternumber;
    end
end
% end NBCoupling
if isfield(GLMType, 'Saccades')
	SAstart = numParams + 1;  SAend = numParams + GLMPars.saccadefilter.filternumber;
	paramind.SA = [SAstart  : SAend];
	numParams = numParams + GLMPars.saccadefilter.filternumber;
end

if GLMType.Contrast
    Cstart = numParams + 1;  Cend = numParams + GLMPars.spikefilters.C.filternumber;
	paramind.C = [Cstart  : Cend];
	numParams = numParams + GLMPars.spikefilters.C.filternumber;
end

% if GLMType.Subunits
%     paramind.SU = (numParams+1):(numParams+GLMPars.subunit_size^2);
%     numParams = numParams+GLMPars.subunit_size^2;
% end


% Assign indices to the stim filter
% If CONVEX assign here
if GLMType.CONVEX
    if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')
        Xstart     = numParams + 1;  
        Xend       = numParams +  GLMPars.stimfilter.frames;
        paramind.X = [Xstart:Xend];
        numParams  = numParams + GLMPars.stimfilter.frames; 
        if isfield(GLMPars.stimfilter,'frames_negative')
            Xstart2     = numParams + 1;  
            Xend        = numParams +  GLMPars.stimfilter.frames_negative;
            paramind.X  = [Xstart:Xend];
            numParams   = numParams + GLMPars.stimfilter.frames_negative;
        end
    elseif strcmp(GLMType.stimfilter_mode, 'fixedSP-ConductanceBased') 
        paramind.Xnote1 = 'Two Filters but fixed spatial filters, excitatory (STA,time1) and inhibitory (STA,time2)';
        Xstart_1 = numParams + 1;  
        Xend_1   = numParams + GLMPars.stimfilter.frames;
        paramind.time1  = [Xstart_1:Xend_1];
        
        Xstart_2 = Xend_1 + 1;  
        Xend_2   = Xend_1 + GLMPars.stimfilter.frames;  
        paramind.time2 = [Xstart_2:Xend_2];
        numParams = numParams  + 2*GLMPars.stimfilter.frames;

        paramind.excitatoryfilter_index = paramind.time1;
        paramind.inhibitoryfilter_index = paramind.time2;
        
        paramind.X = union(paramind.time1,paramind.time2);
    elseif isfield(GLMType,'timefilter') && strcmp(GLMType.timefilter, 'prefilter')
        Xstart = numParams + 1;  
        Xend   = numParams + (GLMPars.stimfilter.ROI_length^2);        
        paramind.X      = [Xstart:Xend];
        paramind.space1 = [Xstart: ((Xstart-1) + (GLMPars.stimfilter.ROI_length^2))];
        numParams       = numParams  +  GLMPars.stimfilter.ROI_length^2; 
    else
        error('you need to properly specifiy the stimfilter in prep_paramind')
    end
end

% If NOT Convex .. Assign here 
if ~GLMType.CONVEX
    convParams = numParams;
    paramind.Xnote0 = 'The stimulus filter X, has non-linear components, breaking convexity';
    
    if strcmp(GLMType.stimfilter_mode, 'rk1') 
        paramind.Xnote1 = 'Stim filter is outer product of space1 and time1';
        Xstart = convParams + 1;  
        Xend   = convParams + (GLMPars.stimfilter.ROI_length^2) + GLMPars.stimfilter.frames;        
        paramind.X      = [Xstart:Xend];
        paramind.space1 = [Xstart: ((Xstart-1) + (GLMPars.stimfilter.ROI_length^2))];
        paramind.time1  = [(Xstart + GLMPars.stimfilter.ROI_length^2) : Xend ];
        numParams       = convParams  +  GLMPars.stimfilter.ROI_length^2 + GLMPars.stimfilter.frames;
    end
    if strcmp(GLMType.stimfilter_mode, 'rk2')
        paramind.Xnote1 = 'Stim filter is sum of outer products of (space1,time1) and (space2,time2)';
        Xstart_1 = convParams + 1;  
        Xend_1   = convParams + (GLMPars.stimfilter.ROI_length^2) + GLMPars.stimfilter.frames;        
        paramind.space1 = [Xstart_1: ((Xstart_1-1) + (GLMPars.stimfilter.ROI_length^2))];
        paramind.time1  = [(Xstart_1 + GLMPars.stimfilter.ROI_length^2) : Xend_1 ];
        Xstart_2 = Xend_1 + 1;  
        Xend_2   = Xend_1 + (GLMPars.stimfilter.ROI_length^2) + GLMPars.stimfilter.frames;  
        paramind.space2 = [Xstart_2: ((Xstart_2-1) + (GLMPars.stimfilter.ROI_length^2))];
        paramind.time2  = [(Xstart_2 + GLMPars.stimfilter.ROI_length^2) : Xend_2 ];
        numParams = convParams  + 2*( GLMPars.stimfilter.ROI_length^2 + GLMPars.stimfilter.frames);
    end
    if  strcmp(GLMType.stimfilter_mode, 'rk2-ConductanceBased')
        paramind.Xnote1 = 'Two Filters, excitatory (space1,time1) and inhibitory (space2,time2)';
        Xstart_1 = convParams + 1;  
        Xend_1   = convParams + (GLMPars.stimfilter.ROI_length^2) + GLMPars.stimfilter.frames;        
        paramind.space1 = [Xstart_1: ((Xstart_1-1) + (GLMPars.stimfilter.ROI_length^2))];
        paramind.time1  = [(Xstart_1 + GLMPars.stimfilter.ROI_length^2) : Xend_1 ];
        Xstart_2 = Xend_1 + 1;  
        Xend_2   = Xend_1 + (GLMPars.stimfilter.ROI_length^2) + GLMPars.stimfilter.frames;  
        paramind.space2 = [Xstart_2: ((Xstart_2-1) + (GLMPars.stimfilter.ROI_length^2))];
        paramind.time2  = [(Xstart_2 + GLMPars.stimfilter.ROI_length^2) : Xend_2 ];
        numParams = convParams  + 2*( GLMPars.stimfilter.ROI_length^2 + GLMPars.stimfilter.frames);

        paramind.excitatoryfilter_index = union(paramind.space1,paramind.time1);
        paramind.inhibitoryfilter_index = union(paramind.space2,paramind.time2);
    end
    
    if GLMType.Subunits
        paramind.SU = (convParams + 1):(convParams + 9);
        % numParams = numParams + 9;
    end
    
    paramind.convParams = convParams;
    paramind.convParams_ind = 1:convParams;
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
   
    