% Version 21
% restarting 2014-03-23
% POST - 
% Lets just do a fixed SP followed by SP filter refit
% NO SPACE in version 21.
% GLMPars.skipframes_eachblock = 120;
% keep it clean .. low on the parameter count

% Version 20
% Add the fixed spatial filter!

% Version 19
%... better graphing
%... better metrics   (do NSEM only for now!) 

% Version 18  DONE 2013-01-16
% Works for NSEM and BW ..  a little hacky on the movie parameters
% not clean yet .. but certainly good enough to work! 
% no good measurements yet.. need to check  for version 19
% Version 18b Starting 2013-11-30
% figure out how to handle NSEM stuff again

% Version 18a  Starting 2013-11-27
% will attempt to end coupling.. more flexible on the 
% will run out of /netapp/snle/lab/temp_matlabcodeAH/glm_AH_18
% Keep all of the directory flaws (mapEI rather than full mapping) 

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



%% Function Calls


% Directories_Params_v20
% postspike_coupling_filterparamsAH
% loadmoviematfile
% findcelltype
% auto_neighbor
% find_paramindGLM
% hack_spatialfilter

%%% Serious downstream changes with new filters %%%
% filterstimulus_train2AH
% linearfilt_grad3AH
%%   Generic initialization 

clear; close all; clear all; clc
% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS

%%%%%% DEFAULT VALUES AS DETERMINED 2013-01-07 !!!!
GLMPars.cone_model = '8pix_Identity_8pix';    GLMPars.cone_sname='p8IDp8';%
%GLMPars.cone_model = '8pix_Model1_1e4_8pix'; GLMPars.cone_sname = 'p8Mod1Max1e4p8';
%GLMPars.cone_model = '8pix_Model1_1e5_8pix'; GLMPars.cone_sname = 'p8Mod1Max1e5p8';
GLMPars.k_filtermode ='rk2';% GLMPars.fixedSPlength = 13;  GLMPars.fixedSP_nullpoint = 'mean'; 
%GLMPars.k_filtermode = 'OnOff_hardrect_fixedSP_STA';% GLMPars.fixedSPlength = 13;  GLMPars.fixedSP_nullpoint = 'mean'; 
GLMPars.fit_type = 'BW';
GLMPars.debug = true;
GLMPars.tolfun = 2; 

GLMPars.func_sname = 'tempmatlabcode/glm21';
GLMPars.binning= 10;       % roughly half a millisec
GLMPars.maxiter = 50;      % 
GLMPars.tolx   = 9;    
         % default 5
%%% Essentially hard wired parameters %%%
GLMPars.psms = 100 ;      %% post spike filter time length in millisecs
GLMPars.cpms = 100 ;      %% cp spike filter time length in millisecs
GLMPars.spcng_psf = pi/2; %% it could be set as pi, but pi/2 is better for "uniform" sampling.
GLMPars.spcng_cp  = pi/2; %% it could be set as pi, but pi/2 is better for "uniform" sampling.
GLMPars.n_psf = 20;       %% or maybe even 243   .. ps basis numbers
GLMPars.K_slen = 13;
GLMPars.STA_Frames = 30;  
GLMPars.n_cp = 8;
GLMPars.Coupling = false;
GLMPars.BiDirect_CP = false;
if ~GLMPars.Coupling
    GLMPars.n_cp = 0;
    GLMPars.BiDirect_Cp = false;
end

GLMPars.map_type = 'mapPRJ';
GLMPars.fullmfilename =mfilename('fullpath');

i_exp = 1;

%%
for i_exp = 1:4
%% Load DataRun, shead_cellID ..  eventually need to be more systematic

    clear exp_nm shead_cellID
    if i_exp == 1, exp_nm = '2012-08-09-3'; shead_cellID = [1471  3676 5086 5161 841 1426 1772 2101 1276];    end
    if i_exp == 2, exp_nm = '2012-09-27-3'; shead_cellID = [1 31 301 1201 1726 91 1909 2360 6858];          end
    if i_exp == 3, exp_nm = '2013-08-19-6'; shead_cellID = [737 1328 1341 2959  2824 3167 3996 5660 6799 5447]; end
    if i_exp == 4, exp_nm = '2013-10-10-0'; shead_cellID = [346 1233 3137 7036 5042 5418 32 768 2778 4354 5866];end
    %              exp_nm = '2013-10-10-0';  head_cellID = [32 768 2778
    %              4354 5866 346 1233 3137 7036 5042 5418]; 

if strcmp(GLMPars.k_filtermode, 'fixedSP')
    GLMPars.K_slen = 1;
end
if strcmp(GLMPars.k_filtermode , 'OnOff_hardrect_fixedSP')
    GLMPars.K_slen = 2;
end
[StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v20(exp_nm, GLMPars.fit_type, GLMPars.map_type)
if ~exist(DirPars.output_dir,'dir'), mkdir(DirPars.output_dir)
%%%%  load cell ids.
if ~exist('shead_cellID','var')
    if strcmp(GLMPars.fit_type , 'BW')
    head_cellID = [datarun_slv.cell_types{2}.cell_ids datarun_slv.cell_types{1}.cell_ids];% datarun_slv.cell_types{5}.cell_ids];
    elseif strcmp(GLMPars.fit_type , 'NSEM')
         [StimulusPars DirPars datarun_slvBW datarun_mas] = Directories_Params_v20_split(exp_nm, 'BW', map_type)
         head_cellIDBW = [datarun_slvBW.cell_types{1}.cell_ids datarun_slvBW.cell_types{2}.cell_ids];     
         head_cellID = [intersect(datarun_slv.cell_types{2}.cell_ids , head_cellIDBW) , intersect(datarun_slv.cell_types{1}.cell_ids , head_cellIDBW) ]
    end
end
if exist('shead_cellID','var')
    head_cellID = shead_cellID;
end

%%%%  Shorten Block count if using Debug
if GLMPars.debug && isfield(StimulusPars, 'BW')
    StimulusPars.BW.FitBlocks = StimulusPars.BW.FitBlocks(1:2);
end
if GLMPars.debug && isfield(StimulusPars, 'NSEM')
    StimulusPars.NSEM.FitBlocks = StimulusPars.NSEM.FitBlocks(1:2);
end
%%% Assigning a SPars
if strcmp(GLMPars.fit_type,   'BW'), SPars = StimulusPars.BW;  SPars.tstim = SPars.default_tstim;   end
if strcmp(GLMPars.fit_type, 'NSEM'), SPars = StimulusPars.NSEM;  end

clear boolean_debug map_type fit_type shead_cellID

%% Setting  save directories, g_opts, init_pars etc.  needs to be cleaned

%%%%%% directory saving happens here%%%%%%%%%%%
if  GLMPars.Coupling,  Main_SolPars = sprintf('%s_ps%d_%s',GLMPars.k_filtermode,GLMPars.n_psf,'cpON');  end
if ~GLMPars.Coupling,  Main_SolPars = sprintf('%s_ps%d_%s',GLMPars.k_filtermode,GLMPars.n_psf,'cpOFF'); end
if  GLMPars.Coupling && GLMPars.BiDirect_CP, Main_SolPars = sprintf('%s_ps%d_%s',GLMPars.k_filtermode,GLMPars.n_psf,'cpBiD'); end

if isfield(GLMPars, 'rect') && GLMPars.rect
    Main_SolPars = sprintf('%s_%s' , GLMPars.rect_type , Main_SolPars );
end          
Other_SolPars = sprintf( 'bin%d_blk%d_tolfun%d',GLMPars.binning,length(SPars.FitBlocks),GLMPars.tolfun);
if GLMPars.Coupling, Other_SolPars = sprintf('%s_cpbasisno%d_neighborpertype%d', Other_SolPars,GLMPars.n_cp, GLMPars.neighbornumber); end
if GLMPars.debug, Other_SolPars = sprintf('debug_%s', Other_SolPars); end
SolutionSpace = sprintf('%s_Fit/%s/%s/%s', GLMPars.fit_type,Main_SolPars, Other_SolPars,GLMPars.cone_model);
d_save = sprintf('%s/%s', DirPars.output_dir, SolutionSpace); 
display(sprintf('Full Model Parameters are:  %s', SolutionSpace)); 
display(sprintf('Save Directory :  %s', d_save));
if ~exist(d_save), mkdir(d_save); end
 clear Main_SolPars Other_SolPars 

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
   'conemodel', GLMPars.cone_model,...
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
   'tstim' , SPars.tstim,...
   'SolutionSpace',SolutionSpace,...
   'd_save', d_save,...
   'map_type', GLMPars.map_type) ;   % ditto

framenum = length(SPars.fitframes)*length(SPars.FitBlocks); % number of frames to train in the long iteration - the results from here are the ones that will ultimatelty be used
init_pars.fitframes = framenum; clear framenum% number of frames to use in the short iterations init_pars.Mk + init_pars.frames -1 is a power of 2 (for fft)
init_pars.fitframes_note = 'fitframes means maximum number of frames ... old legacy name.. will get rid of later';

% number of frames to use to get an initial search point
init_pars.frame_offset = 0;
init_pars.padval = 0;  % something about padding the convolution term 
init_pars.GLMPars      = GLMPars;
init_pars.exp_nm       = DirPars.exp_nm;
init_pars.FULL_StimulusPars = StimulusPars; 
init_pars.DirPars      = DirPars;
init_pars.padval = SPars.PadVal;
init_pars              = postspike_coupling_filterparamsAH(init_pars);
if init_pars.BiDirect_CP && init_pars.Coupling
        init_pars.cp_oneside_filtnumber =   init_pars.cp_filternumber;
        init_pars.cp_filternumber       = 2*init_pars.cp_filternumber;
end
init_pars.cone_sname = GLMPars.cone_sname;
clear framenum SolutionSpace

%% Load Stimulus   -- insert more frame cutting here!
[blockedmoviecell, inputstats, origmatfile] = loadmoviematfile(exp_nm , GLMPars.fit_type, GLMPars.cone_model,'fitmovie');
init_pars.loadedmoviefile = origmatfile;
clear origmatfile;

height = SPars.height;
width  = SPars.width;
fitblocks    = SPars.FitBlocks;
fitframesperblock = length(SPars.fitframes)
totalframes = length(fitblocks) * ( fitframesperblock);
blockstartframe= min(SPars.fitframes);
novelblocks = SPars.NovelBlocks;
concat_fullfitMovie = uint8(zeros(width, height, totalframes));
for i_blk = fitblocks
    blkind = find(fitblocks == i_blk);
    framenums = ( (blkind -1)*fitframesperblock + 1 ) :  (blkind *fitframesperblock);  
    n_blkind = find(novelblocks == i_blk);
    concat_fullfitMovie(:,:,framenums) = blockedmoviecell{n_blkind}.matrix (:,:, blockstartframe:end);    
end
clear height width blkind i_blk n_blkind 
clear staticblocks novelblocks n_blkind i_blk width height blkind 
clear blockedmoviecell blockstartframe fitblocks fitframesperblock framenums
[testmovie] = loadmoviematfile(exp_nm , GLMPars.fit_type, GLMPars.cone_model,'testmovie');
cid = head_cellID(1);  


%%  NEED TO CHANGE THE TIMING WITH WHICH WE PULL SPIKES TO REFLECT CHANGES IN FIT FRAMES!! !! !! !! !!
%%
for cid =  head_cellID
    tic
    %% FIND Neighbors (call auto_neighbor) setup Parmamind
    % make sure is in  Master andSlave dataruns
    clf
    fprintf('now %d of %d..\n',find(head_cellID == cid),length(head_cellID)) ;  % now running on cell 
    
    init_pars.headid         = cid;  
    init_pars.cid            = cid;
    init_pars.note           = 'cid is headid, just trying to get into new naming scheme'
    [Head_CTYPE , head_CTYPE_sname]  = findcelltype(cid, datarun_mas.cell_types);
    init_pars.savename = sprintf('%s_%d', head_CTYPE_sname , cid);
    init_pars.celltype = Head_CTYPE;
    if strcmp(GLMPars.map_type, 'mapEI')
        init_pars.slave_idx         =  find(datarun_slv.cell_ids_map(:,1) == cid);
    end
    if strcmp(GLMPars.map_type, 'mapPRJ')
        init_pars.slave_idx         =  find(datarun_slv.cell_ids == cid);
    end
    init_pars.master_idx        =  find(datarun_mas.cell_ids          == cid);
    
    clear Head_CTYPE head_CTYPE_sname index
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if GLMPars.Coupling
        [ON_Neighbors , OFF_Neighbors ] = auto_neighbor (init_pars , datarun, GLMPars.neighbornumber);
        init_pars.OFF_Neighbor_master  =  OFF_Neighbors;
        init_pars.ON_Neighbor_master   =  ON_Neighbors;
        master_neighbors =[ON_Neighbors , OFF_Neighbors ];
        %slaveidx = [];
        foundneighbors = [];
        lostneighbors  = [];
        for nid = master_neighbors;
           % slaveidx = (find(datarun_slv.cell_ids_map(:,1) == nid));
            if strcmp(GLMPars.map_type, 'mapEI')
                slvidx         =  find(datarun_slv.cell_ids_map(:,1) == nid);
            end
            if strcmp(GLMPars.map_type, 'mapPRJ')
                slvidx      =  find(datarun_slv.cell_ids == nid);
            end       
            if isempty(slvidx);
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
        init_pars.ON_Neighbor   = intersect(  init_pars.ON_Neighbor_master, foundneighbors)

        clear found lost lostneighbors foundneighbors slaveidx lostcell master_neighbors
        clear OFF_Neighbors ON_Neighbors
    end
    if ~GLMPars.Coupling
        init_pars.cp_Neighbors =[];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIALIZING PARAMS %%%%%%%%%%%%%%
    % MAKE PARAM INDEX  %%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [paramind]          =  find_paramindGLM(init_pars)  ;
    init_pars.paramind  = paramind;        
    init_pars.p0        = .01*rand(paramind.numParams,1);
    %init_pars.p0        = 0*rand(paramind.numParams,1);
    
    clear paramind
    %%% Set up Param index   
    
    %% Load OrganizedSpikes (movie and ROI location), initialize Params
    
    %%%%%% Load Spike Times %%%%%%%%%%%%
    clear organizedspikes
	load(sprintf('%s/BlockSpikesandRasters/organizedspikes_%s.mat',DirPars.output_dir,init_pars.savename));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Find firing rate      %%%%%%%
    %%% MU :CRUDE POISSON APPROX %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    totalspikes         = length( cell2mat ( organizedspikes.block.t_sp(SPars.FitBlocks)  ) );
    totalsecs           = length(SPars.FitBlocks) * SPars.seconds_pernovelblock;
    spikespersec        = totalspikes / totalsecs;
    init_pars.recordedfiringrate_hz = spikespersec;  
    init_pars.p0(init_pars.paramind.MU)  = log (spikespersec);
    
    organizedspikes_head = organizedspikes; clear organizedspikes
    clear totalspikes totalsecs spikespersec
    
    
    
    
    %% Construct ROI, grab the Movie for ROI only
    %%% LOAD FIXED SPATIAL FILTR
    %%%%% WORK HERE !!!!
    % Construct the ROI %
    stafit_centercoord = ( datarun_mas.vision.sta_fits{init_pars.master_idx}.mean );
    stafit_sd          = ( datarun_mas.vision.sta_fits{init_pars.master_idx}.sd);
    slvdim.height      = SPars.height; slvdim.width = SPars.width; 
    [center,sd]        = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, StimulusPars.master, slvdim);
    Klen               = GLMPars.K_slen;
    
    
    if strcmp(init_pars.k_filtermode , 'fixedSP' ) || strcmp(init_pars.k_filtermode, 'OnOff_hardrect_fixedSP_STA')
        Klen         = GLMPars.fixedSPlength;
    end
    
    modvec = -floor(Klen/2) : floor(Klen/2) ;
    xdim = center.x_coord + modvec;
    ydim = center.y_coord + modvec;
    if min(xdim)<1
        xdim = 1:Klen;
    end
    if min(ydim)<1
        ydim = 1:Klen;
    end
    
    if max(xdim) > slvdim.width
        xdim = (slvdim.width - Klen + 1):slvdim.width ;
    end
    if max(ydim) > slvdim.height
        ydim = (slvdim.height - Klen + 1):slvdim.height ;
    end
    init_pars.ROI.xdim = xdim;
    init_pars.ROI.ydim = ydim;
   
    %%%%% LOAD UP THE ROI
    novelROIstim_Concat = double(concat_fullfitMovie(xdim,ydim,:));
    if strcmp(GLMPars.fit_type, 'NSEM') 
        novelROIstim_Concat = novelROIstim_Concat / 255 ; 
    end
    if strcmp(GLMPars.fit_type, 'BW') && strcmp(GLMPars.cone_model, '8pix_Model1_1e4_8pix')
        novelROIstim_Concat = novelROIstim_Concat / 255 ; 
    end
    if strcmp(GLMPars.fit_type, 'NSEM') && strcmp(GLMPars.cone_model, '8pix_PtLog_8pix') 
        novelROIstim_Concat = novelROIstim_Concat + 1  ; 
        novelROIstim_Concat = log(novelROIstim_Concat);
        novelROIstim_Concat = novelROIstim_Concat - min(novelROIstim_Concat(:));
        novelROIstim_Concat = novelROIstim_Concat / max(novelROIstim_Concat(:));
    end
    novelROIstim_Concat = reshape(novelROIstim_Concat, Klen^2, init_pars.fitframes);
    if strcmp(init_pars.k_filtermode, 'rk2')
         init_pars.movieoffset_fromzeroone = inputstats.mu_avgIperpix;
         novelROIstim_Concat = novelROIstim_Concat - inputstats.mu_avgIperpix;  
         
         if isfield(GLMPars, 'rect') && GLMPars.rect
            if strcmp(GLMPars.rect_type , 'rect_clipmean')
                if strcmp(head_CTYPE_sname, 'OFFPar')
                    clipind = find(novelROIstim_Concat > 0)
                elseif strcmp(head_CTYPE_sname, 'ONPar')
                    clipind = find(novelROIstim_Concat < 0 )
                end
                
                novelROIstim_Concat(clipind) = 0;
            end
         end
        display(sprintf('~~~ Minimum stim is %d ---' , min(novelROIstim_Concat(:))));
        display(sprintf('~~~ Maximum stim is %d ---' , max(novelROIstim_Concat(:))));
        finalstim = novelROIstim_Concat;
    end

    
    if strcmp(init_pars.k_filtermode, 'fixedSP')
        spfilter = hack_spatialfilter(exp_nm, init_pars.savename, xdim,ydim);       
      %  figure; imagesc(reshape(spfilter,[length(xdim), length(ydim)])); 
        init_pars.spfilter = spfilter;        
        if strcmp( GLMPars.fixedSP_nullpoint , 'mean')
            %init_pars.movieoffset_fromzeroone = mean(novelROIstim_Concat(:));
            %novelROIstim_Concat = novelROIstim_Concat - mean(novelROIstim_Concat(:));               
            init_pars.movieoffset_fromzeroone = inputstats.mu_avgIperpix;
            novelROIstim_Concat = novelROIstim_Concat - inputstats.mu_avgIperpix;                   
        end         
        if isfield(GLMPars, 'rect') && GLMPars.rect
            if strcmp(GLMPars.rect_type , 'rect_clipmean')
                if strcmp(head_CTYPE_sname, 'OFFPar')
                    clipind = find(novelROIstim_Concat > 0);
                elseif strcmp(head_CTYPE_sname, 'ONPar')
                    clipind = find(novelROIstim_Concat < 0 );
                end
                
                novelROIstim_Concat(clipind) = 0;
            end
        end        
        display(sprintf('~~~ Minimum stim is %d ---' , min(novelROIstim_Concat(:))));
        display(sprintf('~~~ Maximum stim is %d ---' , max(novelROIstim_Concat(:))));
        novelROIstim_Concat = (spfilter') * novelROIstim_Concat;  
        finalstim = novelROIstim_Concat;
    end
    
    if strcmp(init_pars.k_filtermode, 'OnOff_hardrect_fixedSP_STA')
        spfilter = hack_spatialfilter(exp_nm, init_pars.savename, xdim,ydim);  
        clf; plot(spfilter);
        init_pars.spfilter_pos = spfilter;
        init_pars.spfilter_neg = spfilter;
        if strcmp( GLMPars.fixedSP_nullpoint , 'mean')            
            init_pars.movieoffset_fromzeroone = inputstats.mu_avgIperpix;
            novelROIstim_Concat = novelROIstim_Concat - inputstats.mu_avgIperpix; 
            posclip = find(novelROIstim_Concat(:) > 0);
            negclip = find(novelROIstim_Concat(:) < 0 );
        end
        pos_stim = novelROIstim_Concat; pos_stim(negclip) = 0;
        neg_stim = novelROIstim_Concat; neg_stim(posclip) = 0;
        display( sprintf('~~ Range for PosClip is [%d,%d] --', min(pos_stim(:)) , max(pos_stim(:)) ));
        display( sprintf('~~ Range for NegClip is [%d,%d] --', min(neg_stim(:)) , max(neg_stim(:)) ));
        pos_stim = (spfilter') * pos_stim;
        neg_stim = (spfilter') * neg_stim;
        
        finalstim.pos = pos_stim;
        finalstim.neg = neg_stim;
    end
    
    clear xdim ydim stafit_centercoord stafit_sd center sd modvec  Klen fullspfilter posclip negclip  totalframes spfilter 
    clear pos_stim neg_stim novelROIstim_Concat slvdim
    %% LOAD Home (and if CP  Neighbor SPIKES)
    if GLMPars.Coupling 
        if isempty(neighbors_id)
            spikesconcat_Neighbors = [];
            
        end
    else
        neighbors_id = [];
        spikesconcat_Neighbors = [];
    end
    
    t_a = SPars.fitseconds(1);
    dur = init_pars.tstim * (length(SPars.fitframes) );
    for nid = [cid, neighbors_id]
        [Head_CTYPE , CTYPE_sname]  = findcelltype(nid, datarun_mas.cell_types);
        savename = sprintf('%s_%d', CTYPE_sname , nid); clear Head_CTYPE CTYPE_sname
        load(sprintf('%s/BlockSpikesandRasters/organizedspikes_%s.mat',DirPars.output_dir,savename));
                
        T_SP = []; blk_count = 0;
        for k = SPars.FitBlocks
            blk_count = blk_count + 1;
            t_sp_full = organizedspikes.block.t_sp_withinblock{k} ; % unit of time: sec, 0 for the onset of the block
            t_sp      = t_sp_full(find(t_sp_full >  t_a));
            t_sp = t_sp - t_a;
            t_spcontext = t_sp + ( blk_count -1 )*dur;
            T_SP = [T_SP ; t_spcontext];
        end
        if nid == cid
            spikesconcat_Home{1} = T_SP;  %%% for some reason when this gets called by glm_ft , putll spikes it needs to be a cell!!!
        else
            ind_neighb = find(neighbors_id == nid);
            spikesconcat_Neighbors{ind_neighb} = T_SP;
        end
        clear organizedspikes t_sp T_SP t_sp_full t_spcontext
    end   
    clear t_a dur blk_count
    clear c N_blk N_fr_blk n_sti_fr n_sp CTYPE ind_neighb isitONParasol isitOFFParasol
    clear CTYPE_sname k neighbors_id nid 


    %% Evaluate the GLM 
    [Basepars,opt_param] = glm_fitAH_21(init_pars,spikesconcat_Home,spikesconcat_Neighbors,gopts,finalstim); % AH stuck in the Comm folder
    %% Saving and Plotting
    opt_param.k_filtermode = Basepars.k_filtermode;
    opt_param.compdate     =  datestr(clock); 
    opt_param.compfunction = GLMPars.func_sname;
  
    Basepars.compfunction = GLMPars.func_sname;
    Basepars.compdate = opt_param.compdate;
    Basepars.fullmfilename = GLMPars.fullmfilename;
    opt_param.fullmfilename = GLMPars.fullmfilename;
    Basepars.opt_param = opt_param;


    xvalperformance = eval_xvalperformance(Basepars, SPars, organizedspikes_head,testmovie{1}, inputstats.mu_avgIperpix);
    opt_param.xvalperformance = xvalperformance;   
    eval(sprintf('save %s/%s.mat Basepars opt_param',Basepars.d_save, Basepars.savename))%,'Track_Progress')
    pdfname = sprintf('%s_%s_%s_%s.pdf', Basepars.fit_type, GLMPars.cone_sname, Basepars.k_filtermode, Basepars.savename);
    printglmfit_AH4(Basepars,opt_param,xvalperformance, Basepars.d_save,pdfname);
    
    


    
end
end

