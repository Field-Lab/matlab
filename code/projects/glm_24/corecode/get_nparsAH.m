function N = get_nparsAH(basepars,neff)

if (~isfield(basepars,'XsqK'))
    basepars.XsqK = 0;
end

if (~isfield(basepars,'meta_filt_mode') && strcmp(basepars.k_filtermode,'stark1')) || strcmp(basepars.k_meta_filt_mode,'starkn')
   %fprintf('note: need to modify codes if starkn means other than stark1\n')
   N = 1+1;
elseif ( strcmp(basepars.k_meta_filt_mode,'indrkn') || strcmp(basepars.k_meta_filt_mode,'stst') )
   N = 1+basepars.stim_n;
elseif (strcmp(basepars.k_meta_filt_mode,'ir6mc') || strcmp(basepars.k_meta_filt_mode,'stc2'))
   N = 1+basepars.stim_n;
elseif (strcmp(basepars.k_meta_filt_mode,'jit2') || strcmp(basepars.k_meta_filt_mode,'stc2r'))
   N = 1+basepars.stim_n;
elseif ( strcmp(basepars.k_meta_filt_mode,'lfs') ||  strcmp(basepars.k_meta_filt_mode,'stc6rs'))
   N = 1+basepars.stim_n;
elseif strcmp(basepars.k_meta_filt_mode,'lfst0')
   N = 1+basepars.stim_n;
elseif strcmp(basepars.k_meta_filt_mode,'lfst')
   N = 1+1+basepars.n_nl;
else
   switch(basepars.k_filtermode)
      case 'sep_basis'
         N = 1+(basepars.XsqK+1)*2*(basepars.nofilters_kspace+basepars.nofilters_ktime);
      case {'sep_raw','rk2'}
         N = 1+(basepars.XsqK+1)*2*(basepars.k_spacepixels+basepars.k_stimframes);
      case {'nonsep','raw'}
         N = 1 + (basepars.XsqK + 1)*basepars.k_spacepixels*basepars.k_stimframes;
      case 'fixfilt'
         N = 1+1;
      case 'fixrk2'
         N = 1+2;
      case {'sta','gbwrk2'}
         N = 1+1;
      case 'gbwrk2c'
         N = 1+1+30;
      case 'gbwrk2mc'
         N = 1+1+30+30;
      case 'lc0'
         N = 1+1+25;
      case {'lcs0','lcs1'}
         N = 1+1+25*2;
      case 'abp'
         N = 1+1+49*2;
   end
end

N = N + basepars.ps_filternumber+neff*basepars.cp_filternumber;

if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
    N = N + length(basepars.ext_timepts);
end

