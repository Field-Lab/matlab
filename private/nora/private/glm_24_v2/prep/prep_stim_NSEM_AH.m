%%%   11-5
%%%  THIS NOW HAS A RESCALING THAT COMES UP STRAIGHT IN FRONT   TESTRUN IS
%%%  PROPOERLY SCALED

%EDOI based on 2011-12-04-0/prep_testrun_eyemov_ROI.m

% AK HEITMAN 2012-09-13  SLIGHT MODIFICATION AND FULL COMMENTING
% NEED TO RUN PRE_STIM_NSEM0 !!!  TO MAKE SURE WE HAVE MOVIE FILES TO USE!!


clear
Directories_Params;
ctype    = 'ON-Parasol'; 
load_cid = 'auto';  %load_cid = 'manual';

%{
% LOAD DATA,  GET DATARUN SET UP
datarun{1}.names.rrs_params_path = sprintf('%s/%s/%s.params', com_dir,dn_mas,dn_mas);
datarun{1}.names.rrs_sta_path    = sprintf('%s/%s/%s.sta',    com_dir,dn_mas,dn_mas);
datarun{1}.default_sta_fits      = 'vision';
datarun{2}.names.rrs_neurons_path = sprintf('%s/%s/%s.neurons', com_dir,dn_slv,dn_slv); 
datarun{2}.names.map_path         = sprintf('%s/%s/ei_map-from-_%s_%s_%s.txt',com_dir,dn_slv, nm_exp, dn_mas, dn_mas);
datarun{2}.default_sta_fits = 'vision';
opt = struct('verbose',1,'load_params',1,'load_neurons',1);
datarun = load_data(datarun,opt);

if isfield(datarun{2}.names,'map_path')
   datarun=load_map(datarun);
else
   datarun=map_cell_types(datarun,'verbose',true);
end

datarun{2}.block.trig = cell(n_blk,1);
for k = 1:n_rep
   trg_oe = datarun{2}.triggers((k-1)*ntb_oe+1:k*ntb_oe);
   datarun{2}.block.trig{2*k-1} = trg_oe(1:ntb_o);
   datarun{2}.block.trig{2*k}   = trg_oe(ntb_o+1:end);
end
%clear dn_* opt
%}
master_height = StimulusPars.master_height
master_width = StimulusPars.master_width;
n_rep         = StimulusPars.n_rep;
n_blk  = StimulusPars.n_blk;
fr_sec = StimulusPars.fr_sec;
BW_NovelBlocks = StimulusPars.BW_NovelBlocks;
BW_StaticBlocks      = StimulusPars.BW_StaticBlocks;
BW_FitBlocks         = StimulusPars.BW_FitBlocks;
BW_PadVal            = StimulusPars.BW_PadVal;
NSEM_NovelBlocks = StimulusPars.NSEM_NovelBlocks;
NSEM_StaticBlocks      = StimulusPars.NSEM_StaticBlocks;
NSEM_FitBlocks         = StimulusPars.NSEM_FitBlocks;
NSEM_PadVal            = StimulusPars.NSEM_PadVal;
NSEM_Refresh           = StimulusPars.NSEM_Refresh;
NSEM_Width             = StimulusPars.NSEM_Width;
NSEM_Height            = StimulusPars.NSEM_Height;
ntb_o                  = StimulusPars.ntb_o;
ntb_e                  = StimulusPars.ntb_e;
ntb_oe                 = StimulusPars.ntb_oe;
nsec_o                 = StimulusPars.nsec_o;
nsec_e                 = StimulusPars.nsec_e;
BW_seedA               = StimulusPars.BW_seedA;
BW_seedC               = StimulusPars.BW_seedC;
BW_seedM               = StimulusPars.BW_seedM;
BW_seedS               = StimulusPars.BW_seedS;
BW_TestFrameIdx        = StimulusPars.BW_TestFrameIdx;
BW_NovelFrameIdx       = StimulusPars.BW_NovelFrameIdx;

%%%%%%%%%%%%%%%%%%%%%%%
psms       = GLMPars.psms
cpms       = GLMPars.cpms;
spcng_psf  = GLMPars.spcng_psf;
spcng_cp   = GLMPars.spcng_cp;
n_psf      = GLMPars.n_psf;
n_cp       = GLMPars.n_cp;
K_slen     = GLMPars.K_slen;
STA_Frames = GLMPars.STA_Frames;
fit_type   = GLMPars.fit_type;




%%%%%
%% LOAD DESIRED CELL ID'S

offpar_prepdir = sprintf('%s/NSEM/OFF-Parasol/prep_STAandROI', output_dir);
onpar_prepdir  = sprintf( '%s/NSEM/ON-Parasol/prep_STAandROI', output_dir);

if strcmp(load_cid, 'auto')    %% INEFFICIENT BUT ONLY TAKES A MINUTE
    CID = [];  
    if strcmp(ctype, 'ON-Parasol')
        CID = datarun{2}.cell_types{1}.cell_ids;
    end
    if strcmp(ctype, 'OFF-Parasol')
        CID = datarun{2}.cell_types{2}.cell_ids;
    end       
end
if strcmp(load_cid, 'manual');
    CID = [ 1 , 31];
end
%% 
for cid = 1:length(CID)
   cid
   testrun.cell_ids = CID(cid);
   testrun.cell_type = ctype;
   testrun.no_trig_bl = [ntb_o, ntb_e];
   testrun.t_trig = datarun{2}.triggers;
   testrun.cell_idx = find(datarun{2}.cell_ids==testrun.cell_ids);
   testrun.t_sp     = datarun{2}.spikes{testrun.cell_idx};
   testrun.master_idx = find(datarun{1}.cell_ids == testrun.cell_ids);
   % note: to check the raster, at this point, use code/2011-12-04-0/check_raster.m
   %-- choose spikes in each block
   datarun{2}.block.t_frame = cell(n_blk,1);
   testrun.block.t_sp = cell(n_blk,1);
   
   for j = 1:n_blk
      datarun{2}.block.t_frame{j} = t_frame_interp(datarun{2}.block.trig{j});
      t_b = datarun{2}.block.t_frame{j}(1);
      t_e = datarun{2}.block.t_frame{j}(end);
      testrun.block.t_sp{j} = testrun.t_sp( (testrun.t_sp > t_b) & (testrun.t_sp < t_e) );
   end
   testrun.block.t_frame = datarun{2}.block.t_frame;
   % eval(sprintf('save %s/testrun_id%d_%s testrun',d_save,testrun.cell_ids,dn_slv));
   
   %-- load visual stimulus
   %ntb = datarun{2}.no_trig_bl; % triggers per block (here 37)
   ts_vs.frames = [30,60]*120; % # frames per block: 30 sec/block x 120 frames/sec
   ts_vs.mov = cell(n_blk,1);
   ts_vs.width = 80;
   ts_vs.height = 40;
   % note: in this experiment, there are only 20 unique stimulus sets, and from 2nd to 20th
   % are used as a unique long-run training data set.
   %datarun{2}.block.trig = cell(100,1);
   d_mov = sprintf('%s/stim_NSEM/Movie',output_dir);
   c = 0;
   for k = [1,2:2:2*59] % 59 unique blocks for non-repeated trials;
      % note (2012-02-20): due to a bug in stim.lisp, every other non-repeated trial
      % has a unique stimulus.  so, instead of 2:2:2*59, 2:4:180 (namely, 45 unique 
      % blocks)
      if k == 1
         c = c+1;
          % Xr: 3600x80x40, [0,255], double.
             load(sprintf('%s/eyemov%d',d_mov,c),'X_rescaled');
          %   Xr = reshape(Xr,3600,80*40); % now 3600x320
             %%% NEW PART    %%%
              X_rescaled = reshape(X_rescaled , size(X_rescaled,1) , size(X_rescaled,2)*size(X_rescaled,3) );
        %     display('here')

         ts_vs.mov{k} = reshape(X_rescaled',80,40,3600);          
      else
         c = c+1;
         Xr0 = zeros(7200,80*40);
         load(sprintf('%s/eyemov%d',d_mov,c),'X_rescaled'); % Xr: 3600x80x40, [0,255], double.
         Xr0(1:3600,   :) = reshape(X_rescaled,3600,80*40); % now 3600x3200
         c = c+1;
         load(sprintf('%s/eyemov%d',d_mov,c),'X_rescaled'); % Xr: 3600x80x40, [0,255], double.
         Xr0(3601:7200,:) = reshape(X_rescaled,3600,80*40); % now 2*3600 x 3200
         
         %%% NEW PART %%%
         %%%%%%%%%%%%%%%%
         
         ts_vs.mov{k} = reshape(Xr0_rescaled',80,40,2*3600); 
      end
   end
   
   % interim summary:
   % now the relevant variables are datarun adn ts_vs.  datarun{2} is relevant
   % in the following analysis (note: datarun{1} is just master).
   
   %-- STA wrt all pixels
   testrun.block.STA  = cell(n_blk,1);
   testrun.block.nST  = zeros(n_blk,1);
   
   for j = 2:2:118 % exclude the 2nd block due to adaptation
      %datarun{2}.block.t_frame{j} = t_frame_interp(datarun{2}.block.trig{j});
      t_b = datarun{2}.block.t_frame{j}(1);
      t_e = datarun{2}.block.t_frame{j}(end);
      testrun.block.t_sp{j} = testrun.t_sp( (testrun.t_sp > t_b) & (testrun.t_sp < t_e) );
      
      IDX      = nan(size(testrun.block.t_sp{j}));
      itv_bin  = mean(diff(datarun{2}.block.trig{j}))/100;
      %-- compute the closest frame index for each spike
      for k = 1:length(testrun.block.t_sp{j})
         t = testrun.block.t_sp{j}(k);
         [min_t,idx] = min(abs(datarun{2}.block.t_frame{j} - t));
         if min_t < itv_bin/2
            IDX(k) = idx;
         end
      end
      IDX = IDX(~isnan(IDX));
      
      testrun.nkt = 30; % number in "k" (stimulus filter) in time (i.e., frames)
      testrun.block.STA{j} = zeros(ts_vs.width,ts_vs.height,testrun.nkt);
      for k = 1:length(IDX)
         idx = IDX(k);
         if idx >= testrun.nkt
            testrun.block.STA{j} = testrun.block.STA{j} + ts_vs.mov{j}(:,:,idx:-1:idx-testrun.nkt+1);
         end
      end
      testrun.block.nST(j) = sum(IDX>=testrun.nkt);
   end
   
   testrun.block.t_frame = datarun{2}.block.t_frame;
   fprintf('DONE\n')
   
   %-- averaging STA computed for each block
   testrun.STA = zeros(size(testrun.block.STA{6}));
   for k = 22:2:118
      testrun.STA = testrun.STA + testrun.block.STA{k};
   end
   
   testrun.STA = testrun.STA/sum(testrun.block.nST);
   
   
   
   %-- note: got to use the index common with the BW
   %[~,max_ind] = max(abs(testrun.STA(:)-1/2)); % max deviation from the mean
   %[max_x,max_y,max_fr] = ind2sub(size(testrun.STA),max_ind);
  %{
   fn_bw = sprintf('%s/output/%s/BW/%s/prep2/testrun_id%d_15sq',...
       glmhome,nm_exp,ctype,testrun.cell_ids);
   trbw = load(fn_bw);
   
   testrun.ROI.o_s = 15; % off-set from the peak % note: 15 is side length, not off-set! (2012-02-20)
   testrun.ROI.ROI_x = trbw.testrun.ROI.ROI_x;
   testrun.ROI.ROI_y = trbw.testrun.ROI.ROI_y;
   testrun.ROI.max_xyfr = trbw.testrun.ROI.max_xyfr;
   clear trbw
   testrun.ROI.STA = testrun.STA(testrun.ROI.ROI_x,testrun.ROI.ROI_y,:);
  %}
   
   max_ind = ( datarun{1}.vision.sta_fits{testrun.master_idx}.mean );
   max_x   = round( max_ind(1)* (ts_vs.width  /master_width)  );
   max_y   = ts_vs.height - round( max_ind(2)* (ts_vs.height /master_height) );
     %%% just a hack .. don't know why
   [~, max_fr] = max(abs(testrun.STA(max_x,max_y,:) - mean(testrun.STA(:)))  );
   
   if rem(K_slen,2) == 1
      testrun.ROI.o_s = (K_slen-1)/2; % off-set from the peak
   end
   %%% RECENTER. . DEFINE BEST K_slen BY K_slen ROI
   [xdim,ydim,zdim] =size ( testrun.STA );
   xhigh = (round(max_x+testrun.ROI.o_s));
   xlow  = (round(max_x-testrun.ROI.o_s));
   if xhigh > xdim
       xhigh = xdim; xlow  = xdim - K_slen + 1;
   end
   if xlow < 1
       xlow = 1; xhigh = K_slen ;
   end 
   yhigh = (round(max_y+testrun.ROI.o_s));
   ylow  = (round(max_y-testrun.ROI.o_s));
   if yhigh > ydim
       yhigh = ydim; ylow  = ydim - K_slen + 1;
   end
   if ylow < 1
       ylow = 1; yhigh = K_slen ;
   end
   testrun.ROI.ROI_x = (xlow:xhigh);
   testrun.ROI.ROI_y = (ylow:yhigh);
   testrun.ROI.max_xyfr = [max_x,max_y,max_fr];
   testrun.ROI.STA = testrun.STA(testrun.ROI.ROI_x,testrun.ROI.ROI_y,:);
   
   
   
   testrun.ROI.mov = cell(100,1);
   for k = [1,2:2:118]
      testrun.ROI.mov{k} = ts_vs.mov{k}(testrun.ROI.ROI_x,testrun.ROI.ROI_y,:);
   end
   
   d_save = sprintf('%s/stim_NSEM/%s/prep_STAandROI',output_dir,ctype);
   if ~exist(d_save , 'dir')
       mkdir(d_save)
   end
   eval(sprintf('save %s/testrun_id%d_ROI testrun',d_save,testrun.cell_ids))
   
end
