% This function computes the gradient of the input term at each time wrt
% the spatial/temporal filters. The index of the neuron in the param vector
% is specified by offset.

% It returns  a 2*(pars.n+pars.Mk) x pars.maxt matrix
% IMPORTANT:
% All of this should be done at the coarser stimulus temporal resolution
% (use Stimpars.dt!)

%%% ALMOST EXPLICITLY FOR RANK-2 DECOMP OF SPAT-TEMP FILTER 
function L = linearfilt_gradAH(p,offset1,Basepars,Stimpars,Trainpars,neuron_no)

[nstim T] = size(Stimpars.movie_ROI);

switch(Basepars.k_filtermode)
                            case {'sep_raw','rk2'}
                                nspace = Basepars.k_spacepixels;
                                ntime = Basepars.k_stimframes;
                            case 'fixfilt' % added, edoi, 2012-01-04
                                nspace = Basepars.k_spacepixels;
                                ntime = Basepars.k_stimframes;
                            case 'sep_basis'
                                nspace = Basepars.nofilters_kspace;
                                ntime = Basepars.nofilters_ktime;
                            case {'nonsep','raw'}
                                error('Cannot get nspace or ntime in nonsep mode!');
                            end

L = zeros(2*(nspace+ntime),T); % t'th column is the gradient of the input term at time t
offset2 = offset1+nspace+ntime;
if (iscell(Basepars.crop_idx))
    space_idx = Basepars.crop_idx{neuron_no};
else
    space_idx = Basepars.crop_idx(:,neuron_no);
end


% Compute gradient of the input term wrt the spatial filters
% The gradient of the input term at time t wrt spatial filter is simply the temporally filtered (vector-valued) stimulus

% Reconvolve using fastconv
% Get the 2 space filters for this neuron
space_filt1 = p(offset1+1:offset1+nspace); % nspace x 1
space_filt2 = p(offset2+1:offset2+nspace); % nspace x 1
% Get the 2 temporal filters for this neuron
temp_filt1 = p(offset1+nspace+1:offset1+nspace+ntime); % ntime x 1
temp_filt2 =p(offset2+nspace+1:offset2+nspace+ntime); % ntime x 1



% If they're parametrized with a basis multiply it out.
if (strcmp(Basepars.k_filtermode,'sep_basis'))
    space_filt1 = Basepars.kspace_basis*space_filt1; % n x 1
    space_filt2 = Basepars.kspace_basis*space_filt2; % n x 1
    temp_filt1 = Basepars.ktime_basis*temp_filt1; % mk x 1
    temp_filt2 = Basepars.ktime_basis*temp_filt2; % mk x 1
end

% Compute gradient of the input term wrt spatial filters
1;
switch(Basepars.k_filtermode)
    case {'sep_raw','rk2'} % edoi
        L(1:nspace,:) = Stimpars.dt*fastconv(Stimpars.movie_ROI(space_idx,:),temp_filt1',nspace,T);
        L(nspace+ntime+1:2*nspace+ntime,:) = Stimpars.dt*fastconv(Stimpars.movie_ROI(space_idx,:),temp_filt2',nspace,T);
    case 'sep_basis'
        1;
        A = Basepars.kspace_basis'*Stimpars.movie_ROI(space_idx,:); % nspace x T
        L(1:nspace,:) = Stimpars.dt*fastconv(A,temp_filt1',nspace,T);
        L(nspace+ntime+1:2*nspace+ntime,:) = Stimpars.dt*fastconv(A,temp_filt2',nspace,T);
end

% Compute gradient of the input term wrt the temporal filters
xs1 = space_filt1'*Stimpars.movie_ROI(space_idx,:); % 1 x T
xs2 = space_filt2'*Stimpars.movie_ROI(space_idx,:); % 1 x T
switch(Basepars.k_filtermode)
    case {'sep_raw','rk2'}
        % The gradient of the input term wrt the temporal filter is the
        % 'stixelized' spatially filtered stimulus
        L(nspace+1:nspace+ntime,ntime:end) = Stimpars.dt .* stim_design(xs1,ntime);
        L(2*nspace+ntime+1:2*(nspace+ntime),ntime:end) = Stimpars.dt .* stim_design(xs2,ntime);
        1;
        for j=1:ntime-1
            L(nspace+1:nspace+j,j) = Stimpars.dt .* xs1(j:-1:1)';
            L(2*nspace+ntime+1:2*nspace+ntime+j,j) = Stimpars.dt .* xs2(j:-1:1)';
        end
    case 'sep_basis'
        % The gradient of the input term at time t wrt kth temporal filter
        % coefficient is jsut the spatially filtered stimulus that is
        % temporally filtered by the kth basis function
        for j=1:ntime
            L(nspace+j,:) = Stimpars.dt.* fastconv(xs1,Basepars.ktime_basis(:,j)',1,T,Basepars.padval);
            L(2*nspace+ntime+j,:) = Stimpars.dt .* fastconv(xs2,Basepars.ktime_basis(:,j)',1,T,Basepars.padval);
        end
end
1;