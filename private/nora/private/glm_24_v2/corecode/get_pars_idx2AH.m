% Function that gets the relevant indices in the parameter vector for a
% given model component

% Mu then Kfilter then PS filter then CP filters   in that order!!
% FM-DEPENDENT (edoi, 2012-01)

function idx = get_pars_idx2AH(Basepars,comp)

offset = 0;
nNeighbors = length(Basepars.cp_Neighbors);

if (strcmp(comp,'ext'))
    if (isfield(Basepars,'ext_timepts') && ~isempty(Basepars.ext_timepts))
        idx = offset+1:offset+length(Basepars.ext_timepts);
        return;
    else
        idx = [];
        return;
        %error('ext_timepts does not exist');
    end
end

%%% DON'T LEAVE THIS PART OUT !!  %%%
if (~isfield(Basepars,'XsqK'))
    Basepars.XsqK = 0;
end

switch(Basepars.k_filtermode)
      case {'nonsep','raw'}
         kfilter_length = Basepars.k_spacepixels*Basepars.k_stimframes;
      case {'sep_raw','rk2'}
         kfilter_length = 2*(Basepars.k_spacepixels+Basepars.k_stimframes);
      case 'sep_basis'
         kfilter_length = 2*(Basepars.nofilters_kspace+Basepars.nofilters_ktime);
      case 'fixfilt'
         kfilter_length = 1;
      case {'sta','gbwrk2'} %,'stark1'
         kfilter_length = 1;
end

stimfilters = kfilter_length( 1 + Basepars.XsqK);

switch (comp)
    case 'mu'
        idx = 1;    
    case 'kfilter'       
        idx = 2:(1+kfilter_length);
    case 'ksq'
        if (~Basepars.XsqK)
            idx = [];
            return;
            %error('XsqK is off or does not exist.');
        end
        
        idx = (1+kfilter_length+1):(1+2*kfilter_length);
    case 'ps'
        idx = (1+stimfilters)+1:(1+stimfilters)+Basepars.ps_filternumber;
    case 'cp'
        idx = (1+stimfilters+Basepars.ps_filternumber) +1 : (1+stimfilters+Basepars.ps_filternumber) +nNeighbors* Basepars.cp_filternumber;        
    case 'kspace'
        [nspace ntime] = get_nspacetime(Basepars);
        idx = [2:1+nspace 1+nspace+ntime+1:1+2*nspace+ntime];
    case 'ktime'
        [nspace ntime] = get_nspacetime(Basepars);
        idx = [1+nspace+1:1+nspace+ntime 1+2*nspace+ntime+1:1+2*(nspace+ntime)];
end

if (isfield(Basepars,'ext_timepts') && ~isempty(Basepars.ext_timepts))
    offset = offset + length(Basepars.ext_timepts);
end
idx = idx + offset;
end












%{
%npars = length(Basepars.p0);
%offset = (j-1)*npars;
if (~isfield(Basepars,'k_meta_filt_mode') && strcmp(Basepars.k_filtermode,'stark1')) || strcmp(Basepars.k_meta_filt_mode,'starkn')
   nkidx = 1;
elseif (strcmp(Basepars.k_meta_filt_mode,'stst') || strcmp(Basepars.k_meta_filt_mode,'indrkn') )
   nkidx = Basepars.stim_n;

elseif (strcmp(basepars.k_meta_filt_mode,'ir6mc') || strcmp(basepars.k_meta_filt_mode,'stc2'))
   nkidx = basepars.stim_n;
elseif (strcmp(basepars.k_meta_filt_mode,'jit2') || strcmp(basepars.k_meta_filt_mode,'stc2r'))
   nkidx = basepars.stim_n;
elseif strcmp(basepars.k_meta_filt_mode,'lfs')
   nkidx = basepars.stim_n;
elseif ( strcmp(basepars.k_meta_filt_mode,'gbwrk2c') || strcmp(basepars.k_meta_filt_mode,'gbwrk2mc'))
   nkidx = basepars.stim_n;   
elseif ( strcmp(basepars.k_meta_filt_mode,'lc0') || strcmp(basepars.k_meta_filt_mode,'lcs0') )
   nkidx = basepars.stim_n;   
elseif ( strcmp(basepars.k_meta_filt_mode,'lcs1') || strcmp(basepars.k_meta_filt_mode,'abp') ) 
   nkidx = basepars.stim_n;
elseif strcmp(basepars.k_meta_filt_mode,'stc6rs')
   nkidx = basepars.stim_n;
elseif strcmp(basepars.k_meta_filt_mode,'lfst0')
   nkidx = basepars.stim_n;
elseif strcmp(basepars.k_meta_filt_mode,'lfst')
   nkidx = basepars.stim_n;
   fprintf('nkidx = basepars.stim_n\n')
%}
        