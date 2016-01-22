%%% Should be the same for each "experiment" .

%%%% Need GLMPars, fit_type, debug_mode, k_filtertype (same as
%%%% GlmPars.k_filtermode)

clear, clear -global, close all

% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS
runtype = 'local';
local_homedir ='/Users/akheitman/Matlab_code/glm_AH_current';
network_homedir ='/netapp/snle/home/snl-e/matlab-standard/code/projects/glm_AH_15'; 
if strcmp(runtype,   'local'), homedir =   local_homedir; end
if strcmp(runtype, 'network'), homedir = network_homedir; end
cd(homedir); path(homedir, path);

%%%%%% DEFAULT VALUES  D NOT CHANGE !!!!
default.GLMPars.binning= 20;       % roughly half a millisec
default.GLMPars.maxiter = 1000;      % 
default.GLMPars.tolx   = 9;    
default.GLMPars.tolfun = 5;          % default 5
default.GLMPars.k_filtermode = 'STA'; 
default.GLMPars.fit_type = 'BW';
default.GLMPars.debug    = false;
default.GLMPars.psms = 100 ;    %% post spike filter time length in millisecs
default.GLMPars.cpms = 100 ;    %%  cp  spike filter time length in millisecs
default.GLMPars.spcng_psf = pi/2; % it could be set as pi, but pi/2 is better for "uniform" sampling.
default.GLMPars.spcng_cp  = pi/2; % it could be set as pi, but pi/2 is better for "uniform" sampling.
default.GLMPars.n_psf = 36;     %%% or maybe even 243   .. ps basis numbers
default.GLMPars.K_slen = 15;
default.GLMPars.STA_Frames = 30;  
default.GLMPars.n_cp = 8;
default.GLMPars.Coupling = true;
default.GLMPars.BiDirect_CP = false;


%%%%%%%%%%%  Standard Important Settings that need to be set
exp_nm                = '2012-08-09-3';  head_cellID = [1471 1786]; %[1276 1426]; %% these will be the lead two for not debugging[1471 1786]
GLMPars               = default.GLMPars;
GLMPars.fit_type      = 'BW';
GLMPars.debug         = true;
GLMPars.k_filtermode  = 'STA';





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Here do your "experiment" where we switch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  variables around 
clear default
for i = 1:3
    if i  == 1
        GLMPars.tolfun = 2;
    elseif i ==2
        GLMPars.tolfun = 3;
    end
    glm_AH_15
end






