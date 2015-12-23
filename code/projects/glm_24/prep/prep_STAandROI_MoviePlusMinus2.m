% 2013-01-9 REVAMP ... Complete OverHaul

% Calls STA_blockcomp .. 
% need to work on the plotting of these guys


% 2012-11-1 UPDATE!!!  BIG CHANGE from .02 .98  to -.5 to .5
 %move the movie back to .5 -.5 in the prep and throughout  


% COMPUTE THE STA, FIND ITS SVD AND INITIALIZE

% 0) INITILIZE  LOAD CELL IDS AND SPIKES
% 1) LOAD THE MOVIE (EITHER NSEM OR BW)
% 2) CYCLE THROUGH CELLS  COMPUTE STA THEN THE ROI
% 3) SAVE INTO TESTRUN STRUCT FOR EACH CELL



%%%%%%%%   0. INITIALIZE: LOAD CELL IDS AND SPIKES     %%%%%%%%%%
%clear
%%
exp_nm = '2012-08-09-3'; boolean_debug = false; map_type = 'mapEI'; 

[StimulusPars DirPars datarun] = Directories_Params_v18(string_date, boolean_debug, 'BW', map_type)   %false means not debug mode




resfactor = 1;
frameshift = 1;
K_slen     = 15;
STA_Frames = 30;




if StimulusPars.BW.commonmovie
    BWmovie_dir = '/snle/home/snl-e/glm/output_AH/BW_Movie';
else
    BWmovie_dir = sprintf('%s/BW/Movie',DirPars.output_dir);  
    if ~exist(BWmovie_dir,'dir'),mkdir(BWmovie_dir),end
end

fit_type = 'BW';
if strcmp(fit_type, 'BW')
    datarun{1} = datarun_all.master;
    datarun{2} = datarun_all.BW_slv;
    SPars = StimulusPars.BW;    
end
if strcmp(fit_type, 'NSEM')
    datarun{1} = datarun_all.master;
    datarun{2} = datarun_all.NSEM_slv;
    SPars = StimulusPars.NSEM;
end

cell_type = 'OFF-Parasol';
CTYPE = {'ON-Parasol','OFF-Parasol'};

%for count = 1:2
%%
output_dir = DirPars.output_dir
for count = 1:2
    if count == 1
        cell_type = 'OFF-Parasol';
    elseif count == 2
        cell_type = 'ON-Parasol';
    end
    d_save = sprintf('%s/%s/%s/prep_STAandROI_res%d_shift%d',output_dir,fit_type,cell_type, resfactor, frameshift);
    if ~exist(d_save,'dir'),mkdir(d_save),end
 %   cd(d_save);
    BW_save = sprintf('%s/%s/%s/prep_STAandROI_2',output_dir,'BW',cell_type);
    NSEM_save=sprintf('%s/%s/%s/prep_STAandROI_2',output_dir,'NSEM',cell_type);
    load_cid = 'auto';    %  load_cid = 'manual';   


    if strcmp(load_cid, 'auto')    %% INEFFICIENT BUT ONLY TAKES A MINUTE
        CID = [];

        if strcmp(cell_type, 'ON-Parasol')
            CID = datarun{2}.cell_types{1}.cell_ids;
        end
        if strcmp(cell_type, 'OFF-Parasol')
            CID = datarun{2}.cell_types{2}.cell_ids;
        end    

    end


%%


    %%%%%%%%%%   1. LOAD MOVIE PARAMS AND MOVIE   use Load Movie  %%%%%%%%%%%%
    %% Movie Loading
    %%% INITIALIZE TS_VS STRUCTURE WHICH WILL HOLD THE MOVIES %% 
    master_height = StimulusPars.master_height;
    master_width = StimulusPars.master_width;
    
    n_rep  = SPars.n_rep;
    n_blk  = SPars.n_blk;
    fr_sec = SPars.fr_sec;
    frames_pertrigger = SPars.frames_pertrigger;
    NovelBlocks = SPars.NovelBlocks;
    FitBlocks = SPars.FitBlocks;
    PadVal    = SPars.PadVal;
    tstim     = SPars.tstim;
    ntb_o = SPars.ntb_o;
    ntb_e = SPars.ntb_e;
    nsec_o = SPars.nsec_o;
    nsec_e = SPars.nsec_e;
    
    if strcmp(fit_type, 'BW')
         TestFrameIdx = SPars.TestFrameIdx;
        NovelFrameIdx = SPars.NovelFrameIdx;
                seedA = SPars.seedA;
                seedC = SPars.seedC;
                seedM = SPars.seedM;
                seedS = SPars.seedS;
    end    




    %%
    %%% LOAD MOVIES INTO TS_VS.M %%
    ts_vs.frames = [nsec_o,nsec_e] * fr_sec; % #FRAMES PER BLOCK .. SET FOR BW STIM
    ts_vs.mov = cell(n_blk,1);  
    switch fit_type
        case 'BW'
           % if ~StimulusPars.BW.commonmovie
           % cd(BWmovie_dir);
            %%% MOVIE SEED PARAMS
            %% for BW fits 60 blocks of 30 seconds .. this takes 1 minute
            %% on alligaor.. un optimized.. 
            sd.a = seedA;    sd.c = seedC;
            sd.m = seedM;    sd.s = seedS; 
            %%% LOAD MOVIES AND RESHAPE MOVIES, SAVE INTO TS_VS STRUCTURE
            %%% THIS BLOCK TAKES SOME TIME .. COUPLE MINUTES
            for k = 1:n_blk     
                %%% KEEP TRACK OF XML FILE NAME, PUT INTO THE DATARUN 
                if rem(k,2) == 1   %%% odd blocks   should be the static movies / test movies 
                    datarun{2}.block.xml{k} = sprintf('%s/xml_interwoven/static/BW-8-1-0.48-11111.xml',BWmovie_dir);
                else
                    sd.s = mod( (sd.a*sd.s + sd.c), sd.m);
                    datarun{2}.block.xml{k} = sprintf('%s/xml_interwoven/novel/BW-8-1-0.48-%s.xml',...
                        BWmovie_dir,int2str(uint32(sd.s)));
                end    
                %%% LOAD THE MOVIE BLOCK INTO MVI STRUCT
                if k == 1
                    mvi = load_movieAH(datarun{2}.block.xml{k}, datarun{2}.block.trig{k});
                    ts_vs.refresh_time = mvi.getRefreshTime;
                    ts_vs.height = mvi.getHeight; 
                    ts_vs.width  = mvi.getWidth;
                    fr_idx = TestFrameIdx ;
                    mov = zeros(ts_vs.width, ts_vs.height, ts_vs.frames(fr_idx));
                    for i = 1:ts_vs.frames(fr_idx)
                            F = mvi.getFrame(i).getBuffer;
                            F = reshape(F,3,ts_vs.width,ts_vs.height);   
                            mov(:,:,i) = F(1,:,:); 
                    end
                    %%%%  BIG CHANGE HERE. CHANGE MOVIE TO PLUS MINUS
                    %%%%  FROM .02 .98  to -.5 +.5
                    mov = mov     - .5;
                    mov = (.5/.48) *mov; 
                    %%%%%%%%%%%%%%%%%%%%%
                    ts_vs.mov{k} = mov; 
                    testmovie    = mov;
                 end
                 if (~isempty( find(NovelBlocks == k)) )
                        mvi = load_movieAH(datarun{2}.block.xml{k}, datarun{2}.block.trig{k});
                        display(sprintf('LOADING BW FIT MOVIE %d out of %d',ceil(k/2),ceil(n_blk/2)) );                 
                        fr_idx = NovelFrameIdx;
                        mov = zeros(ts_vs.width, ts_vs.height, ts_vs.frames(fr_idx));
                        for i = 1:ts_vs.frames(fr_idx)
                            F = mvi.getFrame(i).getBuffer;
                            F = reshape(F,3,ts_vs.width,ts_vs.height);   
                            mov(:,:,i) = F(1,:,:); 
                        end
                        %%%%  BIG CHANGE HERE. CHANGE MOVIE TO PLUS MINUS .5
                        mov = mov     - .5;
                        mov = (.5/.48) *mov; 
                        %%%%%%%%%%%%%%%%%%%%%
                        ts_vs.mov{k} = mov;
                 end
                 if k == 1  && ~StimulusPars.BW.commonmovie
                     display('saving Test Movie')
                     save(sprintf('%s/testmovie.mat', BWmovie_dir),'testmovie'); clear testmovie
                     movieparams = ts_vs;
                     movieparams = rmfield(movieparams, 'mov');
                     save(sprintf('%s/movieparams.mat', BWmovie_dir), 'movieparams');
                 end
                 if (~isempty( find(NovelBlocks == k)) && ~StimulusPars.BW.commonmovie)
                     eval(sprintf('save novelfitmovie_block_%d mov',k));
                 end
                 clear mov; clear mvi;
            end

        
            
            
        case 'NSEM'  %% %  
            %%%%%%%%  STILL NEED TO WORK THIS OUT %%%%%%%%%
            ts_vs.width         = NSEM_Width;
            ts_vs.height        = NSEM_Height;
            ts_vs.refresh_time  = NSEM_Refresh; 
            d_mov = sprintf('%s/NSEM',output_dir);
            NovelBlock = 0;
            for k = [1,NSEM_Blocks] %%% LOAD ALL MOVIES FROM STIM_NSEM FOLDER
                if k == 1
                     load(sprintf('%s/eyemov%d',d_mov,k),'Xr'); % Xr: 3600x80x40, [0,255], double.
                     Xr = reshape(Xr,fr_sec * nsec_o ,NSEM_Width*NSEM_Height); % now 3600x3200
                     ts_vs.mov{k} = reshape(Xr',NSEM_Width,NSEM_Height,fr_sec * nsec_o );  
                     testmovie    = ts_vs.mov{k};
                end

                %%%%%%%%%% THIS SECTION IS A BIT HARD CODED DUE TO HOW THE
                %%%%%%%%%% NSEM MOVIES ARE LOADED AND SAVED INTO MATLAB!!!!
                if ~isempty ( NSEM_Blocks == k )  %%% right now using nec_e is 2x as long as nsec_o
                     display(sprintf('LOADING NSEM FIT MOVIE %d out of %d',ceil(k/2),ceil(n_blk/2)) );   
                     NovelBlock = NovelBlock + 1;
                     Xr0 = zeros(fr_sec*nsec_e,NSEM_Width*NSEM_Height);
                     load(sprintf('%s/eyemov%d',d_mov,NovelBlock),'Xr'); % Xr: 3600x80x40, [0,255], double.
                     Xr0(1:fr_sec * nsec_o,   :) = reshape(Xr,fr_sec*nsec_o,NSEM_Width*NSEM_Height); % now 3600x3200
                     NovelBlock = NovelBlock + 1;
                     load(sprintf('%s/eyemov%d',d_mov,NovelBlock),'Xr'); % Xr: 3600x80x40, [0,255], double.
                     Xr0(fr_sec * nsec_o + 1:fr_sec*nsec_e,:) = reshape(Xr,fr_sec*nsec_o,NSEM_Width*NSEM_Height); % now 2*3600 x 3200
                     ts_vs.mov{k} = reshape(Xr0',NSEM_Width,NSEM_Height,fr_sec*nsec_e); 
                end
            end

            NSEM_dir = sprintf('%s/NSEM',output_dir);
            cd(NSEM_dir);
            save('testmovie.mat','testmovie'); 
            save('tsvs.mat','ts_vs');
            cd(d_save);
    end
    
    
    
    %cd(d_save);




    %% RUN THROUGH EACH CELL .. PER CELL LOOP
    
    for cellnum = 1:length(CID)
    %for cellnum = 1:1
    %% setup   
       clear testrun
       testrun.cell_type = cell_type;
       testrun.cell_ids  = CID(cellnum);
       testrun.no_trig_bl = [ntb_o, ntb_e];
       testrun.t_trig = datarun{2}.triggers;
       testrun.cell_idx = find(datarun{2}.cell_ids==testrun.cell_ids);
       testrun.master_idx = find(datarun{1}.cell_ids == testrun.cell_ids);
       testrun.t_sp     = datarun{2}.spikes{testrun.cell_idx};
       testrun.block.t_sp = cell(n_blk,1);
       testrun.block.t_frame = datarun{2}.block.t_frame;
       


      
      frameshift = 1; Movie = ts_vs.mov; spiketimes = testrun.t_sp;  staframes = 30;
      STA = STA_blockcomp(Movie, staframes, spiketimes, FitBlocks, datarun, resfactor, frameshift);
      figure; imagesc( reshape(STA, [size(STA,1)*size(STA,2), staframes]) ); colorbar
      testrun.STA = STA;
      testrun.frameshift = frameshift;
      testrun.resfactor  = resfactor;
       %%%%%%%%%%%%%%%




       %%%%%%%%
       %% ROI CONSTRUCTION AND STA OVER THE ROI
       %%%%%%%%

       %%% we could also consider the loading center from the master

       %%% (X,Y) CENTER OF MASS  AND THEN ROI OF SIZE 15 15
       %switch fit_type
          % case 'BW'
             % max_ind      = round( datarun{1}.vision.sta_fits{testrun.master_idx}.mean ) % max deviation from the mean   the index of it
              %[t,max_ind2] = max(abs(testrun.STA(:)-1/2)); 
              %[max_x,max_y,max_fr] = ind2sub(size(testrun.STA),max_ind2);
              %maxind2 = [max_x,max_y]

       %%%       
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
       %%% COMPUTE STA OVER JUST THE ROI
       testrun.ROI.mov = cell(n_blk,1);
       for k = [1,FitBlocks]
          if (k == 1 || rem(k,2) == 0)
             testrun.ROI.mov{k} = ts_vs.mov{k}(testrun.ROI.ROI_x,testrun.ROI.ROI_y,:);
          end
       end
       testrun.t_bin = ts_vs.refresh_time/1000; % in sec
       testrun.block.t_frame = datarun{2}.block.t_frame;
       %-- first several significant components in STA (in ROI)
       zmSTA = reshape(testrun.ROI.STA,K_slen^2,STA_Frames);
         % subtract (true) mean value.
       rk = 15;
       [U,S,V] = svds(zmSTA,rk);

       testrun.ROI.zmSTA.STA = zmSTA;
       testrun.ROI.zmSTA.U = U;
       testrun.ROI.zmSTA.S = S;
       testrun.ROI.zmSTA.V = V;

       eval(sprintf('save %s/testrun_id%d_%dsq testrun',d_save,testrun.cell_ids,testrun.ROI.o_s*2+1));
    %%
       %%%%%%%%%%%%%%%%%%
       %%%%% PLOTTING STA MAIN OUTPUT FIGURES
       %%%%%%%%%%%%%%%%%%
       rk=15;
       for k = 1:3  %WE ONLY REALLY NEED THE FIRST 5 PRINCIPAL COMPONENTS 
          figure(1) 
           %%% pause   if you want to watch things change one by one
           clf
          %-- raw STA
          subplot(3,4,1)
          imagesc(zmSTA')
          title(sprintf('STA (ID %4d)',testrun.cell_ids))
          xlabel('Space [stixel]')
          ylabel('Time [frame]')
          colorbar horizontal

          m_fr = testrun.ROI.max_xyfr(3);
          ax = 5:5:15;

          subplot(3,4,3)
          imagesc(reshape(zmSTA(:,m_fr),K_slen,K_slen)), axis image
          set(gca,'xtick',ax), set(gca,'ytick',ax)
          colorbar horizontal
          title(sprintf('Spatial STA at the peak (%d-th frame)',m_fr))

          %-- k-th rank
          idxk = zeros(rk,1);
          idxk(k) = 1;
          STAk = U*diag(idxk.*diag(S))*V';

          subplot(3,4,5)
          imagesc(STAk')
          title(sprintf('%d-th rank component',k))
          xlabel('Space [stixel]')
          ylabel('Time [frame]')
          colorbar horizontal

          m_fr = testrun.ROI.max_xyfr(3);

          subplot(3,4,7)
          imagesc(reshape(STAk(:,m_fr),K_slen,K_slen)), axis image
          set(gca,'xtick',ax), set(gca,'ytick',ax)
          colorbar horizontal
          title(sprintf('Spatial STA (%d-th frame)',m_fr))

          %-- cumulative
          idxc = zeros(rk,1);
          idxc(1:k) = 1;
          STAc = U*diag(idxc.*diag(S))*V';

          subplot(3,4,9)
          imagesc(STAc')
          title(sprintf('Rank-%d STA',k))
          xlabel('Space [stixel]')
          ylabel('Time [frame]')
          colorbar horizontal

          DIFc = zmSTA-STAc;
          subplot(3,4,10)
          imagesc(DIFc')
          err = var(DIFc(:))/var(zmSTA(:))*100;
          title(sprintf('Difference (%2.1f%s)',err,'%'))
          xlabel('Space [stixel]')
          ylabel('Time [frame]')
          colorbar horizontal

          subplot(3,4,11)
          imagesc(reshape(STAc(:,m_fr),K_slen,K_slen)), axis image
          set(gca,'xtick',ax), set(gca,'ytick',ax)
          colorbar horizontal
          title(sprintf('Spatial STA (%d-th frame)',m_fr))

          subplot(3,4,12)
          imagesc(reshape(DIFc(:,m_fr),K_slen,K_slen)), axis image
          xlabel('Space [stixel]')
          ylabel('Time [frame]')
          colorbar horizontal
          title('Difference')

          fl_print = 1;
          if fl_print
             fn = sprintf('STA_cid%d',testrun.cell_ids);
             orient landscape
             if k == 1
                eval(sprintf('print -dpsc2 %s/%s',d_save,fn))
             else
                eval(sprintf('print -dpsc2 -append %s/%s',d_save, fn))
             end
          else
             fprintf('check the figuer (not printing).\n')
             pause
          end
       end
       

    end
end
   
