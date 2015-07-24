function basepars2 = set_nframes(basepars,nframes)

basepars2 = basepars;
basepars2.wholesetsize = nframes; % number of frames to train in the long iteration - the results from here are the ones that will ultimatelty be used
basepars2.maxt = nframes; % number of frames to use in the short iterations init_pars.Mk + init_pars.frames -1 is a power of 2 (for fft)
basepars2.prelimFrames = nframes;  % number of frames to use to get an initial search point

