% Version 17  Capable of an NSEM fit
%               
% Version 16  Solution Space fixed as of 2013-01-07
%             A Bit more efficient
%             Only 1 single spike metric computation - Centi Sec level
%             Move function calling back into glm_AH_16

% Version 15  Network quality code ! ! !
%             some final cleaning of output directory
%             wrapper of , metric stuff as well
%             
% Version 14  my new code on logrpob and spik metrics nice and modular..
%              better version of Directories_params_func_all
%              ready to put on server!!
% Version 13   spikemetrucs code up,  increasingly modular,  better version
% Version 12  .. Automated Neighbor finder  (function auto_neighbor) .. get ready to put onto the
% server .. Bi-Drectional Couplng .. and spike metrics.. cleaner track
% progress .. this is a complete piece of code

% Version 11... Smoothed out the entire process.  Rasters are automatically
%      plotted. FixPS and Fix STA back in nicely.Rasterglm auto processed  

% Version 10 produced some stunning rasters ... atleat for cell1
% Version 10   move the movie back to .5 -.5 in the prep and throughout  
   %  USE OUTPUT_AH_STIMPLUSMINUS
 
% Version 9    serious look at spatial temporal filter.. get rid of Get
       % clean up STA to fixed
       % BIG CHANGE IN INITIALIZATION !!! STIMPARS.DT FACTOR HEREI
% Version 8 fixed spatial temporal filter.. fixed PS filter
% Version 7 is devoted to "FIxing" problems in linear_filt2AH and Add Cross Hess Terms
% Versino 6 is devoted to coupling!!
%without coupling Version 5 worked as of 2012-10-10
% HEAD PROGRAM WHICH LOADS DATA (FROM PREP RUNS) AND EXECUTES THE GLM
% STARTED WITH EDOI@SALK.EDU, 2012-01-16, 2012-04-16.
% AH   begins  June 21  2012  
% THIS ITERATION BEGINS 08-17-2012
% This VERsstion started on   2012-10-3


%%   Generic initialization 
%%%2012-08-09-3    [1471 1786 3676 5086 5161];% OFF
%2012-09-27-3    [1 31 301 1201 1726]; OFF
%%%2012-08-09-3    [841 1426 1772 2101 1276];% ON
%%%2012-09-27-3    [91 1909 2360 6858]; ON
clear; close all
% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS
runtype = 'network';
local_homedir ='/Users/akheitman/Matlab_code/glm_AH_current';
%network_currentdir = '/netapp/snle/home/snl-e/matlab-standard/code/projects/call_glm_AH';
network_homedir ='/netapp/snle/home/snl-e/matlab-standard/code/projects/glm_AH_17_c'; 
if strcmp(runtype,   'local'), homedir =   local_homedir; end
if strcmp(runtype, 'network'), homedir = network_homedir; end
path(homedir, path);
if strcmp(runtype,'network'), path(network_homedir,path); end
cd(homedir); 

%%%%%% DEFAULT VALUES AS DETERMINED 2013-01-07 !!!!
GLMPars.binning= 10;       % roughly half a millisec
GLMPars.maxiter = 5000;      % 
GLMPars.tolx   = 9;    
GLMPars.tolfun = 5;          % default 5
GLMPars.k_filtermode = 'rk2'; 
GLMPars.fit_type = 'NSEM';
GLMPars.debug    = false;
GLMPars.psms = 100 ;    %% post spike filter time length in millisecs
GLMPars.cpms = 100 ;    %%  cp  spike filter time length in millisecs
GLMPars.spcng_psf = pi/2; % it could be set as pi, but pi/2 is better for "uniform" sampling.
GLMPars.spcng_cp  = pi/2; % it could be set as pi, but pi/2 is better for "uniform" sampling.
GLMPars.n_psf = 20;     %%% or maybe even 243   .. ps basis numbers
GLMPars.K_slen = 15;
GLMPars.STA_Frames = 30;  
GLMPars.n_cp = 8;
GLMPars.Coupling = true;
GLMPars.BiDirect_CP = false;
GLMPars.NeighborNumber = 1;
clear local_homedir network_homedir network_currentdir runtype


%%  Left over free variables to play with 
exp_nm = '2012-08-09-3';
head_cellID =  [1471 841 3676 1426 2101 1276 5086 5161] 
% [841 1426 2101 1276] 
%head_cellID = [1471  3676 1786 5086 5161];% OFF
%exp_nm = '2012-09-27-3'; 
%head_cellID = [1 31 301 1201 1726]; % OFF

GLMPars.NeighborNumber = 1;
GLMPars.k_filtermode = 'rk2';
GLMPars.debug = false;
GLMPars.tolfun = 4;
GLMPars.computespikemetrics.BW         = false;
GLMPars.computespikemetrics.NSEM       = false;
GLMPars.computespikemetrics.Error_Only = false;
GLMPars.n_cp = 2;

GLMPars.specialrec = false;


%%
[StimulusPars, DirPars, datarun_all] = Directories_Params_func_all(exp_nm, GLMPars.debug);
DirPars.homedir  = homedir;

if strcmp(GLMPars.fit_type,   'BW'), Slv_StimPars = StimulusPars.BW;   fitdir =   'BW_Fit'; end
if strcmp(GLMPars.fit_type, 'NSEM'), Slv_StimPars = StimulusPars.NSEM; fitdir = 'NSEM_Fit'; end


%%%%%% directory saving happens here%%%%%%%%%%%
if  GLMPars.Coupling,  Main_SolPars = sprintf('%s_ps%d_%s',GLMPars.k_filtermode,GLMPars.n_psf,'cpON');  end
if ~GLMPars.Coupling,  Main_SolPars = sprintf('%s_ps%d_%s',GLMPars.k_filtermode,GLMPars.n_psf,'cpOFF'); end
if  GLMPars.Coupling && GLMPars.BiDirect_CP, Main_SolPars = sprintf('%s_ps%d_%s',GLMPars.k_filtermode,GLMPars.n_psf,'cpBiD'); end

if GLMPars.specialrec, Main_SolPars = sprintf('Rec_%s', Main_SolPars); end
             
Other_SolPars = sprintf( 'ConvTEst_bin%d_blk%d_tolfun%d',GLMPars.binning,length(Slv_StimPars.FitBlocks),GLMPars.tolfun);
if GLMPars.Coupling, Other_SolPars = sprintf('%s_cpneighb_%d', Other_SolPars,GLMPars.NeighborNumber); end
if GLMPars.debug, Other_SolPars = sprintf('debug_%s', Other_SolPars); end
SolutionSpace = sprintf('%s/%s/%s', fitdir,Main_SolPars, Other_SolPars);
d_save = sprintf('%s/%s', DirPars.output_dir, SolutionSpace); 




clear fitdir; clear Main_SolPars Other_SolPars

gopts = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','on',...
   'MaxIter',GLMPars.maxiter,... % you may want to change this
   'TolFun',10^(-(GLMPars.tolfun)),...
   'TolX',10^(-(GLMPars.tolx)));
%%%%% INITIALIZE PARAMETERS %%%%%%
init_pars = struct(...
   'debug_mode', GLMPars.debug,...
   'fit_type',GLMPars.fit_type,...
   'LL_eval_C',false,...
   'maxCGIter',GLMPars.maxiter,...  % CG: conjugate gradient  
   'spikebins_perstimframe',GLMPars.binning,...        % samples/stim frame for post-spike filter (i.e., resolution relative to stim filter)
   'Coupling', GLMPars.Coupling,...
   'BiDirect_CP',GLMPars.BiDirect_CP,...
   'ps_fratio',.4,... %   ???????Pillow's Condition????????  not about resolution but rather about filter length
   'ps_ms',GLMPars.psms,...  % MILLISECS LENGTH OF PS_FILTER  %%% 
   'Nstep', @exp, 'Nprime', @exp, 'Ndoubleprime', @exp,...
   'ps_filternumber',GLMPars.n_psf,...
   'cp_filternumber',GLMPars.n_cp,... % this script is for a single neuron; no coupling of multiple neurons
   'ps_spacing',GLMPars.spcng_psf,...
   'cp_spacing',GLMPars.spcng_cp,...
   'ps_bstretch',0.05,... % time scaling for the raised Gaussian bump
   'ps_alpha',0,...    % onset of ps filter
   'k_filtermode',GLMPars.k_filtermode,...
   'k_meta_filt_mode','rkn',...
   'hessmode','mult',... % 'mult'= quasi-Neuton. this enables a large-scale optimization.
   'cp_fratio',.4,... % ??  PILLLOW CONDITON.. ??
   'cp_ms', GLMPars.cpms,...   %Millisecs of the Couplng Filter
   'cp_alpha',0,...
   'cp_bstretch',0.05,...% time scaling for the raised Gaussian bump
   'k_stimframes' ,   GLMPars.STA_Frames,...
   'old_Gradient',false,...
   'old_Hessian',false,...
   'ps_FIX',  false,...
   'k_spacepixels'  , GLMPars.K_slen^2,...
   'STAscale', 4, ...
   'tstim' , Slv_StimPars.tstim,...
   'SolutionSpace',SolutionSpace,...
   'specialrec', GLMPars.specialrec,...
   'd_save', d_save) ;   % ditto
framenum=Slv_StimPars.fr_sec*Slv_StimPars.nsec_e*length(Slv_StimPars.FitBlocks); % number of frames to train in the long iteration - the results from here are the ones that will ultimatelty be used
init_pars.maxt = framenum; clear framenum% number of frames to use in the short iterations init_pars.Mk + init_pars.frames -1 is a power of 2 (for fft)
% number of frames to use to get an initial search point
init_pars.frame_offset = 0;
init_pars.GLMPars      = GLMPars;
init_pars.exp_nm       = DirPars.exp_nm;
init_pars.FULL_StimulusPars = StimulusPars; 
init_pars.DirPars      = DirPars;
init_pars.padval = Slv_StimPars.PadVal;
init_pars.d_home = homedir;

%%% IMPORTANT  DEFINE BASIS FUNCTINOS OF PS AND COUPLING
init_pars              = postspike_coupling_filterparamsAH(init_pars);
if init_pars.BiDirect_CP && init_pars.Coupling
        init_pars.cp_oneside_filtnumber =   init_pars.cp_filternumber;
        init_pars.cp_filternumber       = 2*init_pars.cp_filternumber;
end
datarun = cell (1,2);
datarun{1} = datarun_all.master;
if strcmp(GLMPars.fit_type, 'BW'),  datarun{2} = datarun_all.BW_slv;  end
if strcmp(GLMPars.fit_type,'NSEM'), datarun{2} = datarun_all.NSEM_slv;end

%%
tic

for cid =  head_cellID
    %% FIND Neighbors (call auto_neighbor) 
    % make sure is in  Master andSlave dataruns
    close all
    index                    = find(head_cellID == cid);
    init_pars.headid         = cid;  
    isitONParasol = ~isempty( find(datarun{1}.cell_types{1}.cell_ids == cid) );
    isitOFFParasol = ~isempty( find(datarun{1}.cell_types{2}.cell_ids == cid) );
    if isitONParasol  && ~isitOFFParasol
        Head_CTYPE = 'ON-Parasol';
    end
    if ~isitONParasol && isitOFFParasol
        Head_CTYPE = 'OFF-Parasol';
    end
    clear isitONParasol isitOFFParasol
    init_pars.celltype = Head_CTYPE;
    init_pars.slave_idx         =  find(datarun{2}.cell_ids_map(:,1) == cid);
    init_pars.master_idx        =  find(datarun{1}.cell_ids          == cid);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       [ON_Neighbors , OFF_Neighbors ] = auto_neighbor2 (cid , datarun{1}, GLMPars.NeighborNumber);
    init_pars.OFF_Neighbor_master  =  OFF_Neighbors;
    init_pars.ON_Neighbor_master   =  ON_Neighbors;
    master_neighbors =[ON_Neighbors , OFF_Neighbors ];
    slaveidx = [];
    foundneighbors = [];
    lostneighbors  = [];
    for nid = master_neighbors;
        slaveidx = (find(datarun{2}.cell_ids_map(:,1) == nid));
        if isempty(slaveidx);
            lostneighbors = [lostneighbors nid];
            display('lostcell');
            lostcell = nid;
        else
            foundneighbors = [foundneighbors nid];
        end
    end
    init_pars.cp_Neighbors  = (foundneighbors); neighbors_id  = foundneighbors;
    init_pars.lostNeighbors = (lostneighbors);
    found                   = foundneighbors;
    lost                    = lostneighbors;   
    init_pars.OFF_Neighbor  = intersect( init_pars.OFF_Neighbor_master, foundneighbors);
    init_pars.ON_Neighbor   = intersect(  init_pars.ON_Neighbor_master, foundneighbors);
    
    
    %%%% hack!!!
    init_pars.ON_Neighbor  = 6946;
    init_pars.OFF_Neighbor = 5986;
    init_pars.cp_Neighbors  = [6946, 5986]
    
    
    clear found lost lostneighbors foundneighbors slaveidx lostcell master_neighbors
    clear OFF_Neighbors ON_Neighbors
    %%%%%%%%%
     
    
    %% Set up Param index   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIALIZING PARAMS %%%%%%%%%%%%%%
    % MAKE PARAM INDEX  %%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(init_pars.k_filtermode , 'raw')
        numParams  = 1 + init_pars.k_stimframes * (init_pars.k_spacepixels) + init_pars.ps_filternumber + length(init_pars.cp_Neighbors)*init_pars.cp_filternumber; 
        MU      = 1;   
        Lstart  = 2;               Lend  = (Lstart-1)  + init_pars.k_stimframes * (init_pars.k_spacepixels);
        PSstart = Lend  + 1;       PSend = (PSstart-1) + init_pars.ps_filternumber;
        CPstart = PSend + 1;       CPend = (CPstart-1) + length(init_pars.cp_Neighbors)*init_pars.cp_filternumber;  
        init_pars.paramind.MU = MU;                MUind  = MU;
        init_pars.paramind.L  = [Lstart  :  Lend]; Lind   = [Lstart  :  Lend];
        init_pars.paramind.PS = [PSstart : PSend]; PSind  = [PSstart : PSend];
        init_pars.paramind.CP = [CPstart : CPend]; CPind  = [CPstart : CPend];   
    elseif strcmp(init_pars.k_filtermode, 'rk2')
        numParams  = 1 + 2*( init_pars.k_stimframes + init_pars.k_spacepixels) + init_pars.ps_filternumber + length(init_pars.cp_Neighbors)*init_pars.cp_filternumber; 
        MU      = 1;   
        Lstart  = 2;                                    Lend  = (Lstart-1)  + 2*( init_pars.k_stimframes + init_pars.k_spacepixels);
        PSstart = Lend  + 1;                            PSend = (PSstart-1) + init_pars.ps_filternumber;
        CPstart = PSend + 1;                            CPend = (CPstart-1) + length(init_pars.cp_Neighbors)*init_pars.cp_filternumber;  
        init_pars.paramind.MU = MU;                     MUind  = MU;
        init_pars.paramind.L  = [Lstart  :  Lend];      Lind   = [Lstart  :  Lend];
        init_pars.paramind.PS = [PSstart : PSend];      PSind  = [PSstart : PSend];
        init_pars.paramind.CP = [CPstart : CPend];      CPind  = [CPstart : CPend]; 
        
        spixels = init_pars.k_spacepixels;     tframes = init_pars.k_stimframes;
        space1_start = Lstart;              space1_end = (space1_start-1) + spixels;
        time1_start  = space1_end + 1 ;      time1_end = (time1_start -1) + tframes;
        space2_start = time1_end+1;         space2_end = (space2_start-1) + spixels;
        time2_start  = space2_end + 1 ;      time2_end = (time2_start -1) + tframes;
        
        init_pars.paramind.SPACE1 = [space1_start : space1_end];  SP1ind   = [space1_start : space1_end];
        init_pars.paramind.TIME1  = [time1_start  :  time1_end];   T1ind   = [time1_start  :  time1_end];
        init_pars.paramind.SPACE2 = [space2_start : space2_end];  SP2ind   = [space2_start : space2_end];
        init_pars.paramind.TIME2  = [time2_start  :  time2_end];   T2ind   = [time2_start  :  time2_end];
	elseif strcmp(init_pars.k_filtermode, 'STA')
        numParams  = 1 + 1 + init_pars.ps_filternumber + length(init_pars.cp_Neighbors)*init_pars.cp_filternumber; 
        MU      = 1;   
        Lscale  = 2;                     
        PSstart = Lscale + 1;                           PSend  = (PSstart-1) + init_pars.ps_filternumber;
        CPstart = PSend + 1;                            CPend  =  (CPstart-1) + length(init_pars.cp_Neighbors)*init_pars.cp_filternumber;  
        init_pars.paramind.MU = MU;                     MUind  = MU;
        init_pars.paramind.L  = Lscale;                 Lind   = Lscale;
        init_pars.paramind.PS = [PSstart : PSend];      PSind  = [PSstart : PSend];
        init_pars.paramind.CP = [CPstart : CPend];      CPind  = [CPstart : CPend]; 
        
        
        if init_pars.specialrec
            init_pars.paramind.SR = CPend + 1;
            numParams = numParams + 1;
        end
        
        
    end
    if ~init_pars.Coupling
        numParams = numParams - length(init_pars.cp_Neighbors)*init_pars.cp_filternumber; 
        init_pars.paramind = rmfield(init_pars.paramind , 'CP' );
    end
    
    %{
    if init_pars.BiDirect_CP && init_pars.Coupling
         CPstart = CPstart;    CPend0 = CPend;
         cpparamsnum0 = 1+ CPend0 - CPstart;
         cpparamsnum  = 2 * cpparamsnum0;
         CPend        = CPstart + cpparamsnum - 1;
         init_pars.paramind.CP = [CPstart : CPend ];
         numParams  =  cpparamsnum0 + numParams;
         CPind = [CPstart: CPend];
    end
    %}
    if init_pars.ps_FIX
        adjust                = init_pars.ps_filternumber - 1;
        init_pars.paramind.PS = init_pars.paramind.PS(1);
        PSind                 = init_pars.paramind.PS(1);
        init_pars.paramind.CP = init_pars.paramind.CP - adjust;
        CPind                 = init_pars.paramind.CP;
        numParams             = numParams - adjust;
    end    
    if init_pars.ps_FIX
        init_pars.ps_Filter  = zeros(init_pars.ps_filternumber,1);
    end
 
    clear  CPend CPstart Lscale MU  PSend  PSstart
    
    %% Load Testrun (movie and ROI location), initialize Params
   % keyboard
    %%%%%% Load BW Testrun %%%%%%%%%%%%
    clear testrun
    fprintf('now %d of %d..\n',index,length(head_cellID))   % now running on cell 
    load(sprintf('%s/%s/%s/prep_STAandROI/testrun_id%d_%dsq',DirPars.output_dir,'BW',Head_CTYPE,cid,GLMPars.K_slen)); 
    
    %%% create large 
    if strcmp(GLMPars.fit_type , 'BW')        
        for k = Slv_StimPars.FitBlocks
            if k == Slv_StimPars.FitBlocks(1)
                novelROIstim_Concat = reshape(testrun.ROI.mov{k}, GLMPars.K_slen^2, (Slv_StimPars.fr_sec * Slv_StimPars.nsec_e ) );
            else
                novelROIstim_Concat = cat(2,novelROIstim_Concat,reshape(testrun.ROI.mov{k},GLMPars.K_slen^2,  (Slv_StimPars.fr_sec * Slv_StimPars.nsec_e ) ));
            end
        end
    end
    
    
    %%% Takes about 2-3 minutes .. relatively cheap
    if strcmp(GLMPars.fit_type , 'NSEM')   
       NSEMpercellROIdir = sprintf('%s/NSEM/ROI_MovieFiles', DirPars.output_dir);
       if (~isdir(NSEMpercellROIdir)),mkdir(NSEMpercellROIdir); end   
       
       if exist(sprintf('%s/%d_ROIMovie.mat',NSEMpercellROIdir, cid))
           display('%%%%%%%%%Loading NSEM video over the ROI from File%%%%%%%%'); 
           load(sprintf('%s/%d_ROIMovie.mat',NSEMpercellROIdir, cid));
       end    
       if ~exist(sprintf('%s/%d_ROIMovie.mat',NSEMpercellROIdir, cid))
           display('%%%%%%%%%Constructing NSEM video over the ROI%%%%%%%%');
           ROI_x = testrun.ROI.ROI_x;
           ROI_y = testrun.ROI.ROI_y;     
           ROIMovie_byblock = cell(1, Slv_StimPars.n_blk);
           
           % LOAD THE STATIC BLOCK
           load (sprintf('/netapp/snle/home/snl-e/glm/output_AH/NSEM_Movie/eyemov%d.mat', 1) ,  'X_rescaled')
           movie_ROI0  = X_rescaled(:,ROI_x,ROI_y);
           movie_ROI2 = reshape( movie_ROI0, [size(movie_ROI0,1), size(movie_ROI0,2)*size(movie_ROI0,3)]);
           movie_ROI = transpose(movie_ROI2);
           
           ROIMovie_byblock{1}.type = 'Static, Raster Generator';
           ROIMovie_byblock{1}.blocknum = 1;
           ROIMovie_byblock{1}.eyemovmatfiles = 1;
           ROIMovie_byblock{1}.movie = movie_ROI;
           ROIMovie_byblock{1}.message = 'All Odd Blocks have the same 30 second Static Movie';
           
           for i_block = Slv_StimPars.NovelBlocks
               %%% 2 files per block during Natural Scenes Fitting
               load (sprintf('/netapp/snle/home/snl-e/glm/output_AH/NSEM_Movie/eyemov%d.mat', i_block ) ,  'X_rescaled')
               movie_ROIa = X_rescaled(:,ROI_x,ROI_y);
               movie_ROIb = reshape( movie_ROIa, [size(movie_ROIa,1), size(movie_ROIa,2)*size(movie_ROIa,3)]);
               movie_ROI1 = transpose(movie_ROIb);  
               clear X_rescaled

               load (sprintf('/netapp/snle/home/snl-e/glm/output_AH/NSEM_Movie/eyemov%d.mat', i_block+1 ) ,  'X_rescaled')
               movie_ROIa = X_rescaled(:,ROI_x,ROI_y);
               movie_ROIb = reshape( movie_ROIa, [size(movie_ROIa,1), size(movie_ROIa,2)*size(movie_ROIa,3)]);
               movie_ROI2 = transpose(movie_ROIb);
               clear X_rescaled
               
               movie_ROI = [movie_ROI1 , movie_ROI2];
               
               ROIMovie_byblock{i_block}.type = 'Novel Fitting';
               ROIMovie_byblock{i_block}.blocknum = i_block;
               ROIMovie_byblock{i_block}.eyemovmatfiles = [i_block, i_block+1];
               ROIMovie_byblock{i_block}.movie = movie_ROI;
               display(sprintf('%d Percent Done', round( (50*i_block) / length(Slv_StimPars.NovelBlocks)) ) );   
               save( sprintf('%s/%d_ROIMovie.mat',NSEMpercellROIdir, cid), 'ROIMovie_byblock');
           end   
       end
       
       %%%    ROIMovie_byblock ...  %%%
       novelROIstim_Concat = [];
       fitblks = Slv_StimPars.FitBlocks;
       for i_block = fitblks
           novelROIstim_Concat = [novelROIstim_Concat , ROIMovie_byblock{i_block}.movie];
       end
    end
     %%%%%%% LOAD MOVIE!! %%%%%%%%%   PROLLy NEEDS SOME fit_type dependence
    % CONCATENATION OF ALL THE FIT BLOCKS. (1-SPACE by TIME)
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% INITIALIZE PARAMETERS %%%%%%%
    %%% MU :CRUDE POISSON APPROX  %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    init_pars.p0        = ones(numParams,1);
    init_pars.p0(1)    = 3;
    %{
 %   init_pars.p0(end) = .1;
    init_pars.p_evolve  = zeros(numParams,init_pars.maxCGIter);
 %   totalspikes         = length( cell2mat
 %   (testrun.block.t_sp(Slv_StimPars.FitBlocks)  ) );  should get changed
 %   
    totalframes         = init_pars.maxt;
    spikespersec        = 100 %(Slv_StimPars.fr_sec * totalspikes ) / totalframes;
    init_pars.p0(MUind) = log (spikespersec);
    PSCP_paramscale     = .2 * log (spikespersec);  % .1 CHOSEN BASED ON PREVIOUS OBSERVATION
    L_paramscale        = init_pars.tstim * log(spikespersec);
    % SET LINEAR FILTER WITH STA, SCALED BY .1*MU
    if strcmp(GLMPars.k_filtermode,'raw')
        K0                 = reshape(testrun.ROI.zmSTA.STA,init_pars.k_spacepixels,init_pars.k_stimframes);
        Lparams0           = reshape(K0 , init_pars.k_spacepixels,init_pars.k_stimframes,1);
        normfactor         = L_paramscale / mean(mean(abs(Lparams0))); 
        init_pars.p0(Lind) = normfactor * Lparams0;
    elseif strcmp(GLMPars.k_filtermode,'rk2')     %%   use same STA from BW as initial guess
        %if strcmp(GLMPars.fit_type , 'BW')
        STA = testrun.ROI.zmSTA.STA;  %NEW AH      i dont' know if this will be useful or not                              
        U   = testrun.ROI.zmSTA.U;  %% zommed in STA SVD   this is only the first 15 dimensions of the STA out of 30
        S   = diag(testrun.ROI.zmSTA.S);
        V   = testrun.ROI.zmSTA.V;          
        %elseif strcmp(GLMPars.fit_type , 'NSEM')
        %    STA = reshape(testrun.ROI.STA, init_pars.k_spacepixels , init_pars.k_stimframes);
        %    [U,S,V]  = svd (STA);
        %    S = diag(S);
        %end
       % S = init_pars.STAscale * S0;
        Lparams0           = [U(:,1)*sqrt(S(1)); V(:,1)*sqrt(S(1)); U(:,2)*sqrt(S(2)); V(:,2)*sqrt(S(2)) ];       
        shift  = min(SP1ind) -1;
        K0                 = Lparams0(SP1ind - shift )*Lparams0(T1ind - shift)' + Lparams0(SP2ind - shift)*Lparams0(T2ind -shift)';
        Lparams0normed     = ( norm(STA)/norm(K0) )  * Lparams0; 
        init_pars.p0(Lind) = init_pars.STAscale *Lparams0normed;
        display('herreinrk2');
        figure; imagesc(STA); colorbar;
        figure; imagesc ( init_pars.STAscale* (norm(STA)/norm(K0)) * K0); colorbar;
        init_pars.STA      = STA;
        
    elseif strcmp(GLMPars.k_filtermode,'STA')
        init_pars.STA      = testrun.ROI.zmSTA.STA;
        init_pars.p0(Lind) = L_paramscale* randn(length(Lind),1);     
    end
    % SET PS AND CP FILTERS.. SCALE NORMAL RV'S
    if ~init_pars.ps_FIX
        init_pars.p0(PSind)     = PSCP_paramscale* randn(length(PSind),1);
    elseif init_pars.ps_FIX
        init_pars.p0(PSind)     = 1;
    end
    if init_pars.Coupling  && ~init_pars.BiDirect_CP
        init_pars.p0(CPind) = PSCP_paramscale*randn(length(neighbors_id)*GLMPars.n_cp,1) ;
    end
    if init_pars.Coupling && init_pars.BiDirect_CP
        init_pars.p0(CPind) = PSCP_paramscale*randn(2*length(neighbors_id)*GLMPars.n_cp,1) ;
    end
    
    
 %   init_pars.p0 = zeros(size(init_pars.p0));
 %   a = init_pars.p0
    %}
    clear totalspikes totalframes spikespersec L_paramscale PSCP_paramscale K0 Lparams0 normfactor
    clear STA U S V shift Lparams0normed
    clear MUind CPind Lind PSind
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % NUMBER OF STIMULUS FRAMES
    % I THINK A CONVERSION TO 8-BIT  FROM .02 TO .98
    %%%% THis probably needs to be writen more cleanly.. less hard coded
%    if strcmp(GLMPars.fit_type ,'NSEM')
%        novelROIstim_Concat = novelROIstim_Concat/255*0.96+0.02; 
%    end


 %%%% Code for adding lock structure for the NSEM  %%%%%
%{
  for j = sta_blocks
                  t_b = datarun{2}.block.t_frame{j}(1);
                  t_e = datarun{2}.block.t_frame{j}(end);
                  % FIND ALL SPIKES WITHIN THE BLOCK
                  testrun.block.t_sp{j} = testrun.t_sp( (testrun.t_sp > t_b) & (testrun.t_sp <= t_e) );             
                  if rem(j,2) == 0 % FOR THE NOVEL FITTING BLOCK
                     IDX      = nan(size(testrun.block.t_sp{j})); %% just initialize vector of size t_sp
                     itv_bin  = mean(diff(datarun{2}.block.trig{j}))/100;  %% really is a constant.. this mean is meaningless haha
                     %-- compute the closest frame index for each spike
                     for k = 1:length(testrun.block.t_sp{j})
                        t = testrun.block.t_sp{j}(k);
                        [min_t,idx] = min(abs(datarun{2}.block.t_frame{j} - t));
                        if min_t < itv_bin/2
                           IDX(k) = idx;
                        end
                     end
                     IDX = IDX(~isnan(IDX));

                     %%% COMPUTE THE STA PER BLOCK
                     testrun.nkt = STA_Frames; % FRAMES IN THE STA  (HOW FAR BACK WE LOOK)
                     testrun.block.STA{j} = zeros(ts_vs.width,ts_vs.height,testrun.nkt);
                     for k = 1:length(IDX)
                        idx = IDX(k);
                        if idx >=testrun.nkt
                           testrun.block.STA{j} = testrun.block.STA{j} + ts_vs.mov{j}(:,:,idx:-1:idx-testrun.nkt+1);
                        end
                     end
                     testrun.block.nST(j) = sum(IDX>=testrun.nkt);  %%  NUMBER OF SPIKES WHICH CONTRIBUTED TO STA CALCULATION
                  end
                  if rem(j,10) == 0
                     fprintf('%3d of %d\r ',j,n_blk)
                  end
        end
%}
    
    %% LOAD NEIGHBOR SPIKES
    if isempty(neighbors_id)
        spikesconcat_Neighbors = [];
    end
    N_blk = length(Slv_StimPars.FitBlocks);
    N_fr_blk = (Slv_StimPars.fr_sec * Slv_StimPars.nsec_e );
    [~, n_sti_fr] = size(novelROIstim_Concat); 
    for nid = [cid, neighbors_id]
        if nid ~= cid
            isitONParasol  = ~isempty( find(datarun{1}.cell_types{1}.cell_ids  == nid) );
            isitOFFParasol = ~isempty( find(datarun{1}.cell_types{2}.cell_ids  == nid) );
            if isitONParasol  && ~isitOFFParasol
                CTYPE = 'ON-Parasol';
            end
            if ~isitONParasol && isitOFFParasol
                CTYPE = 'OFF-Parasol';
            end
        end     
        if nid == cid
            CTYPE = Head_CTYPE;
        end
        if strcmp(GLMPars.fit_type , 'BW')
            load(sprintf('%s/%s/%s/prep_STAandROI/testrun_id%d_%dsq',DirPars.output_dir,GLMPars.fit_type,CTYPE,nid,GLMPars.K_slen)); 
            block = testrun.block;
        elseif strcmp(GLMPars.fit_type , 'NSEM')
           % load(sprintf('%s/%s/%s/prep_STAandROI/testrun_id%d_ROI',DirPars.output_dir,GLMPars.fit_type,CTYPE,nid)); 
            block = datarun{2}.block;
            for i_block = Slv_StimPars.FitBlocks
                
                 slave_idx = find(datarun{2}.cell_ids_map(:,1) == nid);
                 t_sp = datarun{2}.spikes{slave_idx};   
                 t_b = block.t_frame{i_block}(1);
                 t_e = block.t_frame{i_block}(end);
                 block.t_sp{i_block} = t_sp( (t_sp > t_b) & (t_sp <= t_e) );                 
            end
        end
        c = 0;
        n_sp = zeros(2*N_blk,1);
        T_SP = nan(n_sti_fr,1); % pre-define the variable that is large enough
        for k = Slv_StimPars.FitBlocks
            c = c+1;
            t_sp = block.t_sp{k} - block.t_frame{k}(1); % unit of time: sec, 0 for the onset of the block
            n_sp(k) = length(t_sp);
            T_SP(sum(n_sp)-n_sp(k)+1:sum(n_sp)) = t_sp + (c-1)*N_fr_blk*init_pars.tstim; % adding effective onset
        end
        T_SP = T_SP(~isnan(T_SP));
        if nid == cid
            spikesconcat_Home{1} = T_SP;  %%% for some reason when this gets called by glm_ft , putll spikes it needs to be a cell!!!
        else
            ind_neighb = find(neighbors_id == nid);
            spikesconcat_Neighbors{ind_neighb} = T_SP;
        end
        clear testrun t_sp T_SP
    end
    clear c N_blk N_fr_blk n_sti_fr n_sp CTYPE ind_neighb isitONParasol isitOFFParasol
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  INITIALIZE TRACK PROGRESS %%%
    %%%        CALL GLM_FIT        %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Track_Progress = struct('p',cell(GLMPars.maxiter,1),'g',cell(GLMPars.maxiter,1),'H',cell(GLMPars.maxiter,1),'LogProb',0,'counter',0);
    trackdir          = sprintf('%s/Track_Progress',init_pars.d_save);
    if ~exist(trackdir, 'dir'),  mkdir(trackdir); end
    celltrackprog_dir = sprintf('%s/Track_Progress/%d',init_pars.d_save, init_pars.headid);  
    if ~exist(celltrackprog_dir, 'dir'),  mkdir(celltrackprog_dir); end
    save(sprintf('%s/Track_Progress.mat',celltrackprog_dir),'Track_Progress');
    init_pars.trackprog_dir = celltrackprog_dir;
    
    clear celltrackprog_dir trackdir
    %%
    if strcmp(Head_CTYPE, 'ON-Parasol'),  celltype = 'ONPar';  end    
    if strcmp(Head_CTYPE, 'OFF-Parasol'), celltype = 'OFFPar'; end
    fn_save           = sprintf('%s_%d',celltype, cid)
    init_pars.fn_save = fn_save;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% FIT THE MODEL %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 %   cd(DirPars.output_dir);
%    cd(init_pars.SolutionSpace);
    
    [Basepars] = glm_fitAH5_trackprogprob(init_pars,spikesconcat_Home,spikesconcat_Neighbors,gopts,novelROIstim_Concat); % AH stuck in the Comm folder
    printglmfit_AH2_crap(Basepars,Basepars.fn_save, Basepars.d_save );
    computespikemetrics = GLMPars.computespikemetrics;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% 
    simspertrial = 1;
    [~] = glm_FitAccuracy1(Basepars, computespikemetrics, simspertrial);


    
end

toc
