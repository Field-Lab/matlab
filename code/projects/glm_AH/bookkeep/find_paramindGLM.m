% 2014-03-05 akheitman
% hack filled .. bt just make it a function to clear room in the main GLM
% code

function [paramind ] =  find_paramindGLM(Basepars)    


% HACK
if ~(isfield(Basepars, 'cp_Neighbors'))
    Basepars.cp_Neighbors = 0;
end


%%% Set up Param index   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIALIZING PARAMS %%%%%%%%%%%%%%
    % MAKE PARAM INDEX  %%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    paramind.numParams = numParams;
 
    clear  CPend CPstart Lscale MU  PSend  PSstart
    clear space1_end space2_end space1_start space2_start
    clear spixels tframes Lstart Lend SP1ind SP2ind T1ind T2ind
    clear time1_end time1_start time2_start time2_end
    