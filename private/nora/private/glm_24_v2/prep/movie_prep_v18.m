 
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
            for i_blk = 1:n_blk     
 
                %%% LOAD THE MOVIE BLOCK INTO MVI STRUCT
                if i_blk == 1
                    
                    % traditional movie load
                    % movie is naturally loaded as a .02 to .98 file  for noise 
                    tic
                    mvi = load_movieAH(datarun_slv.block.xml{i_blk}, datarun_slv.block.trig{i_blk});
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
                    
                    % consider changing to mov2 = int8(a_mov)
                    toc                    %
                    
                    tic
                    load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/BW-8-1-0.48-11111_RNG_16807/testmovie.mat
                    toc
                    
                    %%%%  BIG CHANGE HERE. CHANGE MOVIE TO PLUS MINUS
                    %%%%  FROM .02 .98  to -.5 +.5
                    mov = mov     - .5;
                    mov = (.5/.48) *mov; 
                    %%%%%%%%%%%%%%%%%%%%%
                    ts_vs.mov{i_blk} = mov; 
                    testmovie    = mov;
                 end
                 if (~isempty( find(NovelBlocks == i_blk)) )
                        mvi = load_movieAH(datarun_slv.block.xml{i_blk}, datarun_slv.block.trig{i_blk});
                        display(sprintf('LOADING BW FIT MOVIE %d out of %d',ceil(i_blk/2),ceil(n_blk/2)) );                 
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
                        ts_vs.mov{i_blk} = mov;
                 end
                 if i_blk == 1  && ~StimulusPars.BW.commonmovie
                     display('saving Test Movie')
                     save(sprintf('%s/testmovie.mat', BWmovie_dir),'testmovie'); clear testmovie
                     movieparams = ts_vs;
                     movieparams = rmfield(movieparams, 'mov');
                     save(sprintf('%s/movieparams.mat', BWmovie_dir), 'movieparams');
                 end
                 if (~isempty( find(NovelBlocks == i_blk)) && ~StimulusPars.BW.commonmovie)
                     eval(sprintf('save novelfitmovie_block_%d mov',i_blk));
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