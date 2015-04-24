% Function that gets the relevant indices in the parameter vector for a
% given model component

% FM-DEPENDENT (edoi, 2012-01)

function idx = get_pars_idxAH(basepars,j,Neff,comp)

npars = get_nparsAH(basepars,Neff);
offset = (j-1)*npars;

if (strcmp(comp,'ext'))
    if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
        idx = offset+1:offset+length(basepars.ext_timepts);
        return;
    else
        idx = [];
        return;
        %error('ext_timepts does not exist');
    end
end

if (~isfield(basepars,'XsqK'))
    basepars.XsqK = 0;
end

if (~isfield(basepars,'k_meta_filt_mode') && strcmp(basepars.k_filtermode,'stark1')) || strcmp(basepars.k_meta_filt_mode,'starkn')
   nkidx = 1;
elseif (strcmp(basepars.k_meta_filt_mode,'stst') || strcmp(basepars.k_meta_filt_mode,'indrkn') )
   nkidx = basepars.stim_n;
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
else
   switch(basepars.k_filtermode)
      case {'nonsep','raw'}
         nkidx = basepars.k_spacepixels*basepars.k_stimframes;
      case {'sep_raw','rk2'}
         nkidx = 2*(basepars.k_spacepixels+basepars.k_stimframes);
      case 'sep_basis'
         nkidx = 2*(basepars.nofilters_kspace+basepars.nofilters_ktime);
      case 'fixfilt'
         nkidx = 1;
%      case {'fixrk2','indrk2'} % edoi
%         nkidx = 2;
      case {'sta','gbwrk2'} %,'stark1'
         nkidx = 1;
   end
end

switch (comp)
    case 'b'
        idx = 1;    
    case 'k'
        
        idx = 2:(1+nkidx);

    case 'ksq'
        if (~basepars.XsqK)
            idx = [];
            return;
            %error('XsqK is off or does not exist.');
        end
        
        idx = (1+nkidx+1):(1+2*nkidx);

    case 'ps'
        
        idx = 1+(basepars.XsqK+1)*nkidx+1:(1+(basepars.XsqK+1)*nkidx+basepars.ps_filternumber);
        
    case 'cp'
        idx = (1+(basepars.XsqK+1)*nkidx+basepars.ps_filternumber+1):(1+(basepars.XsqK+1)*nkidx+basepars.ps_filternumber+(Neff-1)*basepars.cp_filternumber);
        
    case 'kspace'
        [nspace ntime] = get_nspacetime(basepars);
        idx = [2:1+nspace 1+nspace+ntime+1:1+2*nspace+ntime];
    case 'ktime'
        [nspace ntime] = get_nspacetime(basepars);
        idx = [1+nspace+1:1+nspace+ntime 1+2*nspace+ntime+1:1+2*(nspace+ntime)];
end

if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
    offset = offset + length(basepars.ext_timepts);
end

idx = idx + offset;
        