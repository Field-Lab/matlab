lGrad = Trainpars.lGrad;  %%% FROM LINEARFILT2_AH

kx = kx;  %%%FULL FILTERED STIMULUS FROM FILTERSTIMULUS_TRAIN2AH

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


nspace = Basepars.k_spacepixels;
ntime  = Basepars.k_stimframes;

ind1 = 1                      : nspace;
ind2 = (nspace+1)             : (nspace + ntime);
ind3 = (nspace+ntime+1)       : 2*nspace + ntime;
ind4 = (2*nspace + ntime + 1) : 2* (nspace + ntime);


guess2 =( temp_filt1')  *  (lGrad(ind2,:)) + ...
        ( temp_filt2')  *  (lGrad(ind4,:)) ;
     
guess1 = (space_filt1')  *  (lGrad(ind1,:)) + ...
        (space_filt2')  *  (lGrad(ind3,:));

    