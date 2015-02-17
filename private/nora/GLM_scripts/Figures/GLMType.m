function GLMType=GLMType()

%GLMType.stimfilter_mode = 'rk1';
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.CONVEX = true;
GLMType.CouplingFilters = true;
GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
GLMType.map_type = 'mapPRJ';
GLMType.debug = false;
GLMType.specialchange = false;
GLMType.CBP=false;
GLMType.input_pt_nonlinearity      = false;
GLMType.input_pt_nonlinearity_type = 'piece_linear_aboutmean';
GLMType.DoubleOpt = false;
GLMType.TonicDrive = true;
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.Subunits = false;
GLMType.func_sname = 'glmwrap24_CP';
GLMType.fullmfilename =mfilename('fullpath');
GLMType.fit_type='NSEM';

end