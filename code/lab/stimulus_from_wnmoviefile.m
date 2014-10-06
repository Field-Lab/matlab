function stimulus = stimulus_from_wnmoviefile(wnmoviefile)

wnmovie = wnmoviefile.getWhiteMovieParams();
stimulus.field_height = wnmovie.height;
stimulus.field_width  = wnmovie.width;
stimulus.seed         = wnmovie.seed;
stimulus.constrast    = wnmovie.contrastValue;
stimulus.sparse       = wnmovie.sparse;
stimulus.probability  = wnmovie.probability;
stimulus.colortype    = char(wnmovie.colorType);
stimulus.movietype    = char(wnmovie.movieType);
stimulus.rng          = char(wnmovie.rng);


% FIXME: copy paste from STIMULUS_FROM_JAVA_MOVIE should be abstracted
switch stimulus.colortype
    case 'INDEPENDENT'
        stimulus.independent = 't';
    case 'DEPENDENT'
        stimulus.independent = 'nil';
    case 'SEPARATED'
        stimulus.independent = 't';
        stimulus.separated = 't';
    otherwise
        error('Color type ''%s'' not recognized.',movie.getUnderlyingMovie.colorType.toString.toCharArray')
end

% FIXME: copy paste from STIMULUS_FROM_JAVA_MOVIE should be abstracted
switch stimulus.movietype
    case 'GAUSSIAN_MOVIE'
        stimulus.type = ':gaussian';
    case 'BINARY_MOVIE'
        stimulus.type = ':binary';
    otherwise
        error('Noise type ''%s'' not recognized.',movie.getUnderlyingMovie.noiseType.toString.toCharArray')
end