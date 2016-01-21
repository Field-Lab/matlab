% 2014-02-11   adding fixedSp as an option

% JUST GET RID OF THE STIMPARS.DT

% This function computes the gradient of the input term at each time wrt
% the spatial/temporal filters. The index of the neuron in the param vector
% is specified by offset.

% It returns  a 2*(pars.n+pars.Mk) x pars.maxt matrix
% IMPORTANT:

%%% AH: I THINK THE BELOW COMMENT IS WRONG!!!
% All of this should be done at the coarser stimulus temporal resolution
% (use Stimpars.dt!)

%%% ALMOST EXPLICITLY FOR RANK-2 DECOMP OF SPAT-TEMP FILTER 
function L = linearfilt_grad2AH(p,Basepars,Stimpars,Trainpars,neuron_no)

%[nstim T] = size(Stimpars.movie_ROI);
[spacePixels , stimFrames] = size(Stimpars.movie_ROI);


if strcmp(Basepars.k_filtermode , 'fixedSP')
    nspace = 1;
    ntime  = 5%Basepars.k_stimframes
    L = zeros(ntime,stimFrames);
    display('running stim design');
    xs = Stimpars.movie_ROI;
    
    
    L(: , ntime:end) = stim_design(xs,ntime);
	for j=1:ntime-1
        L(1: j , j) =  xs(j:-1:1)';
    end   
end





if ~strcmp(Basepars.k_filtermode , 'fixedSP')
    switch(Basepars.k_filtermode)
                                case {'sep_raw','rk2'}
                                    nspace = Basepars.k_spacepixels;
                                    ntime  = Basepars.k_stimframes;
                                case 'fixfilt' % added, edoi, 2012-01-04
                                    nspace = Basepars.k_spacepixels;
                                    ntime  = Basepars.k_stimframes;
                                case 'sep_basis'
                                    nspace = Basepars.nofilters_kspace;
                                    ntime  = Basepars.nofilters_ktime;
                                case {'nonsep','raw'}
                                    error('Cannot get nspace or ntime in nonsep mode!');
                                end

    L      = zeros(2*(nspace+ntime),stimFrames); % t'th column is the gradient of the input term at time t
    % PARAM INDICES
    s1_idx = Basepars.paramind.SPACE1 ;
    t1_idx = Basepars.paramind.TIME1  ;
    s2_idx = Basepars.paramind.SPACE2 ;
    t2_idx = Basepars.paramind.TIME2  ;

    % Get the 2 space filters for this neuron
    space_filt1 = p(s1_idx); % nspace x 1
    space_filt2 = p(s2_idx); % nspace x 1
    % Get the 2 temporal filters for this neuron
    temp_filt1 = p(t1_idx); % ntime x 1
    temp_filt2 = p(t2_idx); % ntime x 1

    % If they're parametrized with a basis multiply it out.
    if (strcmp(Basepars.k_filtermode,'sep_basis'))
        space_filt1 = Basepars.kspace_basis * space_filt1; % n x 1
        space_filt2 = Basepars.kspace_basis * space_filt2; % n x 1
        temp_filt1  = Basepars.ktime_basis  *  temp_filt1; % mk x 1
        temp_filt2  = Basepars.ktime_basis  *  temp_filt2; % mk x 1
    end

    % Compute gradient of the input term wrt spatial filters
    switch(Basepars.k_filtermode)
        case {'sep_raw','rk2'} % edoi
            L(1:nspace,:)                      = fastconvAH(Stimpars.movie_ROI,temp_filt1',nspace,stimFrames);
            L(nspace+ntime+1:2*nspace+ntime,:) = fastconvAH(Stimpars.movie_ROI,temp_filt2',nspace,stimFrames);
        case 'sep_basis'
            A = Basepars.kspace_basis'*Stimpars.movie_ROI(space_idx,:); % nspace x stimFrames
            L(1:nspace,:)                      = fastconvAH(A,temp_filt1',nspace,stimFrames);
            L(nspace+ntime+1:2*nspace+ntime,:) = fastconvAH(A,temp_filt2',nspace,stimFrames);
    end

    % Compute gradient of the input term wrt the temporal filters
    xs1 = space_filt1'*Stimpars.movie_ROI; % 1 x stimFrames
    xs2 = space_filt2'*Stimpars.movie_ROI; % 1 x stimFrames
    switch(Basepars.k_filtermode)
        case {'sep_raw','rk2'}
            %%% runs it into a  30 x stimFrames matrix
            L( nspace+1         : nspace+ntime    , ntime:end) = stim_design(xs1,ntime);
            L( 2*nspace+ntime+1 : 2*(nspace+ntime), ntime:end) = stim_design(xs2,ntime);
            for j=1:ntime-1
                L( nspace+1         : nspace+j        , j) =  xs1(j:-1:1)';
                L( 2*nspace+ntime+1 : 2*nspace+ntime+j, j) =  xs2(j:-1:1)';
            end
        case 'sep_basis'
            % The gradient of the input term at time t wrt kth temporal filter
            % coefficient is jsut the spatially filtered stimulus that is
            % temporally filtered by the kth basis function
            for j=1:ntime
                L(nspace+j,:)         = fastconvAH(xs1,Basepars.ktime_basis(:,j)',1,stimFrames,Basepars.padval);
                L(2*nspace+ntime+j,:) = fastconvAH(xs2,Basepars.ktime_basis(:,j)',1,stimFrames,Basepars.padval);
            end
    end
    if isfield(Basepars,'old_Gradient')
        if Basepars.old_Gradient 
            L = Stimpars.dt * L;
        end
    end
end
