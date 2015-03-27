% function      movieinfo = NSEM_rawmovieinformation(string_NSEMmoviename)
%               2013-12-14  AKHeitman 
%
% usage:        Holds the base parameters/information of our ever expanding list of NSEM .rawMovies 
%               Information that should be  independent of Lisp call 
%
% arguments:    string_NSEMmoviename
%
%
% calls:        none
%
% outputs:      movieinfo .. basic information independent of the lisp call      
%
% paths:        run glmpath_18.m before everything



function  movieinfo = NSEM_rawmovieinfo(string_NSEMmoviename)
% Lay down the parameters here
rawmovie_name = string_NSEMmoviename;
default_tstim = .00832750;
switch rawmovie_name
    case 'eye-120-3_0-3600'
 
        movieinfo.moviename               = rawmovie_name;
        movieinfo.totalframes             = 120*3600; 
        movieinfo.frames_persecond        = 120;
        movieinfo.totalseconds            = 3600;
        movieinfo.ideal_refreshrate       = default_tstim; 
        movieinfo.frames_perimage         = 120;
        movieinfo.dimension               = [320 160];
        movieinfo.CRTbasepixel            = '2by2';
        movieinfo.ideal_rr_note           = 'should be the measured frame rate in seconds when taking triggers into account';
        movieinfo.creatingmatfile         = 'lost';
        movieinfo.creator                 ='MGreschner';
        movieinfo.create_date             ='late 2011?';
        movieinfo.rawmoviedir_SALK        = 'somewher in rush/snlearchive/ ..stimuli/ ...';
        movieinfo.note                    = 'main go to NSEM stim until mid 2013';
        movieinfo.rawmoviefile            = ''; 
        
    case  'eye-long-v2'
        movieinfo.moviename               = rawmovie_name;
        movieinfo.totalframes             = 120*3600; 
        movieinfo.frames_persecond        = 120;
        movieinfo.totalseconds            = 3600;
        movieinfo.ideal_refreshrate       = default_tstim;
        movieinfo.frames_perimage         = 120;
        movieinfo.dimension               = [320 160];
        movieinfo.CRTbasepixel            = '2by2';
        movieinfo.ideal_rr_note           = 'should be the measured frame rate in seconds when taking triggers into account';
        movieinfo.creatingmatfile         = '?';
        movieinfo.creator                 ='MGreschner and AKHeitman';
        movieinfo.create_date             ='~August 2013';
        movieinfo.rawmoviefile_SALK       = '/braid/snle/data/temporary_stimuli/eye-long-v2.rawMovie';
        movieinfo.rawmoviefile            = '/braid/snle/data/temporary_stimuli/eye-long-v2.rawMovie';
        
    case 'eye-traces'   
        movieinfo.moviename               = rawmovie_name;
        movieinfo.totalframes             = 120*40; 
        movieinfo.frames_persecond        = 120;
        movieinfo.totalseconds            = 40;
        movieinfo.ideal_refreshrate       = default_tstim; 
        movieinfo.frames_perimage         = 120;
        movieinfo.dimension               = [320 160];
        movieinfo.CRTbasepixel            = '2by2';
        movieinfo.creator                 ='MGreschner and AKHeitman';
        movieinfo.ideal_rr_note           ='should be the measured frame rate in seconds when taking triggers into account';
        movieinfo.creatingmatfile         = '?';
        movieinfo.create_date             ='~August 2013?';
        movieinfo.rawmoviefile_SALK       = '/braid/snle/data/temporary_stimuli/eye-traces.rawMovie';
        movieinfo.rawmoviefile            = '/braid/snle/data/temporary_stimuli/eye-traces.rawMovie';
        
	case 'FEM900FF_longrast' 
        movieinfo.moviename               = rawmovie_name;
        movieinfo.totalframes             = 120*(62*60); 
        movieinfo.frames_persecond        = 120;
        movieinfo.totalseconds            = 3600;
        movieinfo.ideal_refreshrate       = default_tstim; 
        movieinfo.frames_perimage         = 120;
        movieinfo.dimension               = [320 160];
        movieinfo.CRTbasepixel            = '2by2';
        movieinfo.creator                 = 'AKHeitman and MGreschner';
        movieinfo.ideal_rr_note           = 'should be the measured frame rate in seconds when taking triggers into account';
        movieinfo.creatingmatfile         = '?';
        movieinfo.create_date             = '~August 2013?';
        movieinfo.rawmoviefile_SALK       = '/Akheitman/VanHateren/newrawmovie/FEM900FF_longrast.rawMovie  on Alligator';
        movieinfo.rawmoviefile            = '/Akheitman/VanHateren/newrawmovie/FEM900FF_longrast.rawMovie';
        movieinfo.note                    = 'in vision this is mislabeled as 439200 frames.. it is really 446400 frames';
        movieinfo.note2                    = 'embeeded a FEM versus no FEM test within the raster... as well as strong stim to weak stim for reconstruction';


end

%{
rawmovie_name = 'eye-traces.rawMovie'; plots = true;
framesperblock = 4800;
reduced_dim = [80 40];
blocks = 1;
params.moviename      = rawmovie_name;
params.reduced_dim    = reduced_dim;
params.framesperblock = framesperblock;
%params.blocks         = blocks;
params.secondsperblock= 100;
params.refreshrate    = '120hz';
params.noteX          = 'unaltered movie.. in uint8 form.. goes from 0 to 256 I belive.. it is a linear scale';
params.noteXr         = 'downsampled movie  to a 80 by 40... equivalent to standard 8 pixel block used in white noise runs';
params.noteX_rescaled = 'everything shifted to a +- pt.5 scale for easier mathematical interpretation';
params.noteX          = 'all movies are index1 = frame, index2 = xval, index3 = yval';
FOI = 120:120:framesperblock;  plotsperblock = 40; picsperplot = 1; subplot1 = 1; subplot2 = 1;
%}

%{
rawmovie_name = 'eye-traces-2.rawMovie'; plots = true;
framesperblock = 12000;
reduced_dim = [80 40];
blocks = 1;
params.moviename      = rawmovie_name;
params.reduced_dim    = reduced_dim;
params.framesperblock = framesperblock;
%params.blocks         = blocks;
params.secondsperblock= 100;
params.refreshrate    = '120hz';
params.noteX          = 'unaltered movie.. in uint8 form.. goes from 0 to 256 I belive.. it is a linear scale';
params.noteXr         = 'downsampled movie  to a 80 by 40... equivalent to standard 8 pixel block used in white noise runs';
params.noteX_rescaled = 'everything shifted to a +- pt.5 scale for easier mathematical interpretation';
params.noteX          = 'all movies are index1 = frame, index2 = xval, index3 = yval';
FOI = 120:120:framesperblock;  plotsperblock = 25; picsperplot = 4; subplot1 = 2; subplot2 = 2;
%}

%{
rawmovie_name = 'eye-traces-2b.rawMovie'
plots = true
framesperblock = 12000;
reduced_dim = [80 40];
blocks = 1;
params.moviename      = rawmovie_name;
params.reduced_dim    = reduced_dim;
params.framesperblock = framesperblock;
params.blocks         = blocks;
params.secondsperblock= 30;
params.refreshrate    = '120hz';
params.noteX          = 'unaltered movie.. in uint8 form.. goes from 0 to 256 I belive.. it is a linear scale';
params.noteXr         = 'downsampled movie  to a 80 by 40... equivalent to standard 8 pixel block used in white noise runs';
params.noteX_rescaled = 'everything shifted to a +- pt.5 scale for easier mathematical interpretation';
params.noteX          = 'all movies are index1 = frame, index2 = xval, index3 = yval';
FOI = 120:120:framesperblock;  plotsperblock = 25; picsperplot = 4; subplot1 = 2; subplot2 = 2;
%}
end       
