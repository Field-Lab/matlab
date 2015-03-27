%%%  AK HEITMAN   2012-12-17
%%%2012-08-09-3    [1471 1786 3676 5086 5161];%
%%%2012-09-27-3    [1 31 301 1201 1726]

%%% Template for future SOl Space "experiments"
%%% Make sure parallel processes never touch the same folder or cells
%%% 3 cells from each data set ! ! 
%%% 
clear; clear -global; close all
runtype = 'local';

% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS
local_homedir ='/Users/akheitman/Matlab_code/glm_AH_current';
network_currentdir = '/netapp/snle/home/snl-e/matlab-standard/code/projects/call_glm_AH';
network_homedir ='/netapp/snle/home/snl-e/matlab-standard/code/projects/glm_AH_15b'; 
if strcmp(runtype,   'local'), homedir =   local_homedir; end
if strcmp(runtype, 'network'), homedir = network_homedir; end
path(homedir, path);
if strcmp(runtype,'network'), path(network_homedir,path); end
cd(homedir); 


%%%%%% DEFAULT VALUES  D NOT CHANGE !!!!
default_GLMPars.binning= 20;       % roughly half a millisec
default_GLMPars.maxiter = 1000;      % 
default_GLMPars.tolx   = 9;    
default_GLMPars.tolfun = 5;          % default 5
default_GLMPars.k_filtermode = 'STA'; 
default_GLMPars.fit_type = 'BW';
default_GLMPars.debug    = false;
default_GLMPars.psms = 100 ;    %% post spike filter time length in millisecs
default_GLMPars.cpms = 100 ;    %%  cp  spike filter time length in millisecs
default_GLMPars.spcng_psf = pi/2; % it could be set as pi, but pi/2 is better for "uniform" sampling.
default_GLMPars.spcng_cp  = pi/2; % it could be set as pi, but pi/2 is better for "uniform" sampling.
default_GLMPars.n_psf = 36;     %%% or maybe even 243   .. ps basis numbers
default_GLMPars.K_slen = 15;
default_GLMPars.STA_Frames = 30;  
default_GLMPars.n_cp = 8;
default_GLMPars.Coupling = true;
default_GLMPars.BiDirect_CP = false;

clear local_homedir network_homedir network_currentdir runtype

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CAN START CHANGING HERE !!! %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% STANDARD IMPORTANT SETTINGS THAT NEED TO BE SET AND SUBJECT TO CHANGE
GLMPars               = default_GLMPars;
GLMPars.fit_type      = 'BW';
GLMPars.debug         = false;
GLMPars.k_filtermode  = 'STA';
exp_nm  = '2012-08-09-3';  head_cellID = [1471];%1786]; % [1276 1426];
GLMPars.tolfun = 2;
GLMPars.computespikemetrics.BW         = true;
GLMPars.computespikemetrics.NSEM       = true;
GLMPars.computespikemetrics.Error_Only = true;

%%% DO YOUR "EXPERIMENT" HERE
tic
for i = 1
    if i  == 1
        GLMPars.tolfun = 2;
    end
    GLMPars
    glm_AH_15
end
toc