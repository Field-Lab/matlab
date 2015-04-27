function [Basepars] = glm_fitAH6_trackprogprob(Basepars,spikesconcat_Home,spikesconcat_Neighbors,gopts,novelROIstim_Concat)
%

% called by glm_AH_5

% This is a rip-off version of chaitu's mutiple_neuron.m.  Here the code is
% simplified to fit parameters of a *single* neuron (i.e., no coupling term
% between neurons).
% edoi@salk.edu, 2012-01-16,2012-04-16.

%calls create_histbasis
%calls pull_spikes
%calls grad_basis
%calls train_glm2 (rk2)
%calls nonsep_lgrad(raw)   train_glm_nonsep('raw'(

% new work and variable renaming for now

% NOT AN ACTUAL FUNCTION .. JUST CLEARS SPACES




%% Setup Trainpars
% LOAD LOGICAL SPIKES INTO INTO TRAINPARS
% CONVOLVE SPIKES WITH BASIS FUNCTIONS
% MAKE ADJUSTMENT IN THE CASE OF BI-DIRECTIONAL COUPLING

Trainpars.dt = Basepars.tstim / Basepars.spikebins_perstimframe; % time resolution of spikes    
Trainpars.neighbors_id = Basepars.cp_Neighbors; Trainpars.baseneuron_idx = 1; 
neighbors_id = Basepars.cp_Neighbors;
microBins = Basepars.spikebins_perstimframe*Basepars.maxt;  % this is in (1/50)th units per frame
duration  = Basepars.tstim * Basepars.maxt;
if (~isfield(Basepars,'frame_offset'))
   stim_offset = 0;
else
   stim_offset = Basepars.tstim*Basepars.frame_offset;
end

%%% PULL SPIKES   %%%%%%%
display('%%% Loading Spikes into Logical on Microbin Scale %%%')
Trainpars.logicalspike_microbin_Home      = sparse(logical(false(microBins,1)));         %in some sense initialize b4 puling spikes
Trainpars.logicalspike_microbin_Neighbors = sparse(logical(false(microBins,length(neighbors_id)))); 
[Trainpars.sptimes_Home     ,Trainpars.logicalspike_microbin_Home     ,Trainpars.negSpikes_Home     ] = pull_spikesAH(spikesconcat_Home,Basepars.headid,microBins,duration,Trainpars.dt,stim_offset);
if Basepars.Coupling && length(Basepars.cp_Neighbors) >0
    [Trainpars.sptimes_Neighbors,Trainpars.logicalspike_microbin_Neighbors,Trainpars.negSpikes_Neighbors] = pull_spikesAH(spikesconcat_Neighbors,neighbors_id,microBins,duration,Trainpars.dt,stim_offset);
end

display('%%% Convolving spikes with PS and CP basis vectors %%%%%%');
microBins_offset = Basepars.frame_offset*Basepars.spikebins_perstimframe;
Trainpars.psbasisGrad = grad_basisAH([Trainpars.negSpikes_Home; Trainpars.logicalspike_microbin_Home],Basepars.ps_basis,microBins_offset);
Trainpars.cpbasisGrad = [];
if Basepars.Coupling && Basepars.BiDirect_CP
    oldNeighspike = Trainpars.logicalspike_microbin_Neighbors;
    couplingoffset = Basepars.cp_timebins;
    newCPpart1    = oldNeighspike( (couplingoffset+1):end , : );
    newCPpart2    = false( couplingoffset , size(newCPpart1,2) );
    newCP         = [newCPpart1 ; newCPpart2 ];
    Trainpars.logicalspike_microbin_Neighbors_preshift = oldNeighspike ; 
    Trainpars.logicalspike_microbin_Neighbors = newCP;   
end
if Basepars.Coupling && length(Basepars.cp_Neighbors >0)
    Trainpars.cpbasisGrad = grad_basisAH([Trainpars.negSpikes_Neighbors; Trainpars.logicalspike_microbin_Neighbors],Basepars.cp_basis,microBins_offset);% note: coupling is not examined here.
end

clear microBins duration stim_offset neighbors_id microBins_offset
clear couplingoffset newCPpart1 newCPpart2 newCP oldNeighspike

%% Setup Stimpars  final changes to Basepars for wierd ooptions
%keyboard
%Nice and concatedneate   stim,   taking account for any offset
Stimpars.dt = Basepars.tstim;
frame_idx = Basepars.frame_offset+1:Basepars.frame_offset+Basepars.maxt;   %%% total frames.. again   
Stimpars.movie_ROI = novelROIstim_Concat(:,frame_idx); clear frame_idx  %%% the movie over the ROI

%%%%%%%%%%%% FIXED STA  or FIXED PS FILTER   %%%%%%%%%%%%%%%%%
if strcmp(Basepars.k_filtermode, 'STA')
    %%% scalar function of time %%%%
    Basepars.kx_STA          = sum(fastconv(Stimpars.movie_ROI,Basepars.STA,Basepars.k_spacepixels,Basepars.maxt,Basepars.padval),1)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Dirty Code for Rectification %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Basepars.specialrec
        onefilter = ones(size(Basepars.STA,1) , 10); 
        mov = Stimpars.movie_ROI;
        %%% m_coit   = movie convolved over the intergration time  ~8 msec
        m_coit =  fastconv(mov,onefilter,Basepars.k_spacepixels,Basepars.maxt,Basepars.padval);
        pos_parts = find(m_coit >= 0) ;
        neg_parts = find(m_coit < 0 ) ;
        if strcmp(Basepars.celltype, 'ON-Parasol'), rec_m_coit = m_coit; rec_m_coit(neg_parts) = 0; end     
        if strcmp(Basepars.celltype, 'OFF-Parasol'), rec_m_coit = m_coit; rec_m_coit(pos_parts) = 0; end

        K = Basepars.STA;
        [~,mx_idx] = max(abs(K(:)));
        [~,~,mx_fr] = ind2sub([15,15,30],mx_idx);
        Spatial_Filter = Basepars.STA(:,mx_fr);
        
        Basepars.rect = rec_m_coit;
        Basepars.rect_conv_spSTA=  sum(fastconv(rec_m_coit, Spatial_Filter, Basepars.k_spacepixels,Basepars.maxt,Basepars.padval),1)';
        figure; subplot(3,1,1); imagesc( reshape(Spatial_Filter, 15,15) ); colorbar;  title('spatial filter from STA peak')
        subplot(3,1,3); hist(Basepars.rect_conv_spSTA); 
        subplot(3,1,2); hist(Basepars.rect(:),50);
        eval(sprintf('print -dpdf %s/Rec_STA_%s.pdf',Basepars.d_save, Basepars.fn_save) )        
end


if Basepars.fullrec
        
        t_filter = ones(size(Basepars.STA,1),1) * Basepars.int_time_filter'; 
        mov = Stimpars.movie_ROI;
        %%% m_coit   = movie convolved over the intergration time
        m_coit =  fastconv(mov,t_filter,Basepars.k_spacepixels,Basepars.maxt,Basepars.padval);
        pos_parts = find(m_coit >= 0) ;
        neg_parts = find(m_coit < 0 ) ;
        if strcmp(Basepars.celltype, 'ON-Parasol'), rec_m_coit = m_coit; rec_m_coit(neg_parts) = 0; end     
        if strcmp(Basepars.celltype, 'OFF-Parasol'), rec_m_coit = m_coit; rec_m_coit(pos_parts) = 0; end
        
        sqrec = rec_m_coit .^2;
        
        if strcmp(Basepars.celltype, 'OFF-Parasol'), sqrec = - sqrec; end
      
        Basepars.rect_full = sqrec;
        figure; subplot(2,1,1); plot(Basepars.int_time_filter); title('integratointime kernel')
        %subplot(3,1,3); hist(Basepars.rect_conv_spSTA); 
        subplot(2,1,2); hist(Basepars.rect_full(:),50); title('hist of rectified integrated movie')
        eval(sprintf('print -dpdf %s/Rec_STA_%s.pdf',Basepars.d_save, Basepars.fn_save) )        
end


    
if Basepars.ps_FIX
    ps_weights               = Basepars.ps_Filter;
    PS                       = Basepars.ps_basis*ps_weights; % Basepars.Mhist x N    
    Basepars.ps_FIXconvolved = (ps_weights') *cell2mat(Trainpars.psbasisGrad);
    clear ps_weights PS
end

%% Actually Run the Optimization and Save
p0 = Basepars.p0;
Basepars.crop_idx = (1:Basepars.k_spacepixels)';
timeoffset = Basepars.frame_offset*Stimpars.dt; % offset time in seconds
display('%%% Starting fminunc %%%%');
switch Basepars.k_filtermode
   case 'rk2'
     % [p_opt f_opt g_opt,H_opt,exitflag,output] = train_glm2_AH(p0,Basepars,Stimpars,Trainpars,gopts,1);   %%% all good till here
       pstar  = p0;
        if isfield(Basepars, 'rect_conv_spSTA')
            [pstar fstar eflag output] = fminunc(@(p) ll_func3AH(p,Basepars,Stimpars,Trainpars,[],Basepars.rect_conv_spSTA),pstar,gopts);
           [f_opt g_opt H_opt lcifs] = ll_func3AH(pstar,Basepars,Stimpars,Trainpars,[], Basepars.rect_conv_spSTA);
       else
            [pstar fstar eflag output] = fminunc(@(p) ll_func3AH(p,Basepars,Stimpars,Trainpars),pstar,gopts);
            [f_opt g_opt H_opt lcifs] = ll_func3AH(pstar,Basepars,Stimpars,Trainpars);
        end
       p_opt = pstar;
    case 'STA'
       pstar  = p0;
       
       if isfield(Basepars, 'rect_conv_spSTA')
           [pstar fstar eflag output] = fminunc(@(p) ll_func3AH(p,Basepars,Stimpars,Trainpars,Basepars.kx_STA,Basepars.rect_conv_spSTA),pstar,gopts);
           [f_opt g_opt H_opt lcifs] = ll_func3AH(pstar,Basepars,Stimpars,Trainpars,Basepars.kx_STA, Basepars.rect_conv_spSTA);
       else
           [pstar fstar eflag output] = fminunc(@(p) ll_func3AH(p,Basepars,Stimpars,Trainpars,Basepars.kx_STA),pstar,gopts);
           [f_opt g_opt H_opt lcifs] = ll_func3AH(pstar,Basepars,Stimpars,Trainpars,Basepars.kx_STA);
       end
       p_opt = pstar;
   case 'raw'  % raw = nonsep; pixel- and frame-based filter, i.e., no assumption on K.
      Trainpars.lgrad = nonsep_lgrad(Basepars,stimpars,Trainpars);
      [p_opt,f_opt,g_opt,H_opt,exitflag,output] = train_glm_nonsep(p0,Basepars,Stimpars,Trainpars,gopts);
end
Basepars.p_opt  = p_opt;
opt_param.p     = p_opt;
opt_param.p0    = p0;
opt_param.f     = f_opt;
opt_param.g     = g_opt;
opt_param.H     = H_opt;
opt_param.exitflag = eflag;
opt_param.output   = output;
opt_param.lcifs = lcifs;
Basepars.p_opt  = p_opt;   %%% Final Change to Basepars... this is what gets saved ! !
%Stimpars = rmfield(Stimpars,'movie_ROI');
Trainpars = rmfield(Trainpars,'psbasisGrad');
if(isfield(Trainpars,'lgrad'))
   Trainpars = rmfield(Trainpars,'lgrad');
end
if(isfield(Trainpars,'cpbasisGrad'))
   Trainpars = rmfield(Trainpars,'cpbasisGrad');
end
cid = Basepars.headid;

%load(sprintf('%s/Track_Progress.mat',Basepars.trackprog_dir));
save(sprintf('%s/%s.mat',Basepars.d_save, Basepars.fn_save), 'Basepars','opt_param','Trainpars','cid','eflag','output')%,'Track_Progress'); 
%print_glm_fitAH(cid,opt_param,Basepars,Trainpars,Basepars.fn_save,Track_Progress);
clear eflag output Track_Progress lcifs p_opt p0 f_opt g_opt H_opt pstar fstar
clear timeoffset

%delete(sprintf('%s/Track_Progress.mat',Basepars.trackprog_dir));

end

