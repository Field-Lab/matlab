%%
homedir ='/Users/akheitman/Matlab_code/glm_AH_current';
path(homedir, path);

clear Track_Progress
Directories_Params; clear datarun; clear dn_slv;
cd(Basepars.d_save);
cd ..
if (~isdir('Rasters'))
    mkdir('Rasters'); 
end
cd Rasters;
stim_fit = Basepars.fit_type;

for count = 1:2
    %% RELOAD DATARUN FOR SAKE OF GRABBING
    clear datarun;  %%% get rid of datarun and start again
    if count ==1 
        stim_vrf = 'BW';
    end
    if count ==2
        stim_vrf = 'NSEM';
    end
    if strcmp(stim_vrf    , 'NSEM')
        dn_slv = dn_slvNSEM;
    elseif strcmp(stim_vrf , 'BW')
        dn_slv = dn_slvBW;
    end    
    %%% TRIGGERS ASSOCIATED WITH STIM TYPE, SINCE STIM TYPE GETS MODULATED %%%
    if strcmp(stim_vrf,'BW');
       StimulusPars.ntb_o = 13;  % TRIGGERS PER STATIC RASTER BLOCK
       StimulusPars.ntb_e = 37;  % TRIGGERS PER NOVEL FITTING BLOCK

       StimulusPars.nsec_o = 10; % SECONDS PER STATIC RASTER BLOCK  
       StimulusPars.nsec_e = 30; % SECONDS PER NOVEL FITTING BLOCK
    end
    if strcmp(stim_vrf,'NSEM')
       StimulusPars.ntb_o = 37;  % TRIGGERS PER STATIC RASTER BLOCK
       StimulusPars.ntb_e = 73;  % TRIGGERS PER NOVEL FITTING BLOCK 

       StimulusPars.nsec_o = 30; % SECONDS PER STATIC RASTER BLOCK  
       StimulusPars.nsec_e = 60; % SECONDS PER NOVEL FITTING BLOCK
    end
    StimulusPars.ntb_oe = StimulusPars.ntb_o+StimulusPars.ntb_e;

    %%% SETUP DATARUN FOR LOADING THE MASTER DATA
    datarun{1}.names.rrs_params_path  = sprintf('%s/%s/%s.params', com_dir,dn_mas,dn_mas);
    datarun{1}.names.rrs_sta_path     = sprintf('%s/%s/%s.sta',    com_dir,dn_mas,dn_mas);
    datarun{1}.default_sta_fits       = 'vision';

    %%% SETUP DATARUN FOR ENSLAVED (BW OR NSEM .. CAN'T CLASSIFY ON OWN DUE TO MOVIE ISSUES)
    datarun{2}.names.rrs_neurons_path = sprintf('%s/%s/%s.neurons', com_dir,dn_slv,dn_slv);
    datarun{2}.names.map_path         = sprintf(... 
       '%s/%s/cellmatch_mapEI_from_%s.txt',com_dir,dn_slv,dn_mas);
    datarun{2}.default_sta_fits = 'vision';

    %%% ACTUALLY LOAD UP DATA
    opt = struct('verbose',1,'load_params',1,'load_neurons',1);
    datarun = load_data(datarun,opt);
    if isfield(datarun{2}.names,'map_path')
        datarun=load_map(datarun);
    else
        datarun=map_cell_types(datarun,'verbose',true);
    end
    clear opt

    %%% ADD BLOCK STRUCTURE TO THE ENSLAVED DATA
    datarun{2}.block.trig = cell(StimulusPars.n_blk,1);    
    for k = 1:StimulusPars.n_rep
       trg_oe = datarun{2}.triggers((k-1)*StimulusPars.ntb_oe+1:k*StimulusPars.ntb_oe);   %  this loop fills in the time of the triggers
       datarun{2}.block.trig{2*k-1} = trg_oe(1:StimulusPars.ntb_o);
       datarun{2}.block.trig{2*k}   = trg_oe(StimulusPars.ntb_o+1:end);
    end
    datarun{2}.block.t_frame = cell(StimulusPars.n_blk,1);
    for j = 1:StimulusPars.n_blk
          datarun{2}.block.t_frame{j} = t_frame_interp(datarun{2}.block.trig{j});
          t_b = datarun{2}.block.t_frame{j}(1);
          t_e = datarun{2}.block.t_frame{j}(end);
    end
    datarun{2}.block.xml  = cell(StimulusPars.n_blk,1);   %% empty cells for now
    datarun{2}.block.trig = cell(StimulusPars.n_blk,1);   %%% empty cells for now
    for k = 1:StimulusPars.n_rep
        trg_oe = datarun{2}.triggers((k-1)*StimulusPars.ntb_oe+1:k*StimulusPars.ntb_oe);
        datarun{2}.block.trig{2*k-1} = trg_oe(1:StimulusPars.ntb_o);
        datarun{2}.block.trig{2*k}   = trg_oe(StimulusPars.ntb_o+1:end);
    end

    %%  LOAD THE MOVIES AND TESTRUN 
    cid        = Basepars.headid;
    master_idx = find(datarun{1}.cell_ids           == cid);
    slave_idx  = find(datarun{2}.cell_ids_map(:,1)  == cid);  
    %%%%  LOADING TESTRUN %%%%%%%%%
    load(sprintf('%s/%s/%s/prep_STAandROI/testrun_id%d_%dsq',output_dir,GLMPars.fit_type,Basepars.celltype,cid,GLMPars.K_slen)); 

        
    %%%% LOADING MOVIES %%%%%%%
    if count == 1
        testmovie_ROI = testrun.ROI.mov{1};
        frames = size(testmovie_ROI, 3);
        testmovie_ROI = reshape(testmovie_ROI, Basepars.k_spacepixels , frames);
    elseif count == 2
        movies_dir = '/snle/home/snl-e/glm/output_AH_STIMplusminus';
        K_slen = GLMPars.K_slen;
        eval(sprintf('load %s/NSEM_Movie/eyemov1.mat',movies_dir));
        testmovie0 = X_rescaled;
        clear X; clear Xr; clear X_rescaled;
        movie_width  = StimulusPars.NSEM_Width;
        movie_height = StimulusPars.NSEM_Height;
        testmovie1 = reshape(testmovie0, size(testmovie0,1), movie_width * movie_height);
        testmovie2 = transpose(testmovie1);
        testmovie  = reshape(testmovie2, movie_width, movie_height, size(testmovie0,1));
        
        ROI_Xind = testrun.ROI.ROI_x;
        ROI_Yind = testrun.ROI.ROI_y;
        testmovie_ROI = testmovie( ROI_Xind , ROI_Yind , : );
        frames = size(testmovie_ROI, 3);
        testmovie_ROI = reshape(testmovie_ROI, Basepars.k_spacepixels , frames);
        clear testmovie0; clear testmovie1; clear testmovie2;
    end
    % put into form which is easy for downstream programs to read
    Stimpars.movie_ROI  = testmovie_ROI;   clear testmovie_ROI; clear testmovie;
    frames = size(Stimpars.movie_ROI , 2);

    %% LOAD SPIKES  FOR THE RETINA RASTERS
    block.t_sp = cell(StimulusPars.n_blk,1);
    t_sp     = datarun{2}.spikes{slave_idx};
    for j = 1:StimulusPars.n_blk
             datarun{2}.block.t_frame{j} = t_frame_interp(datarun{2}.block.trig{j});   %interp stands for interpolation
             t_b = datarun{2}.block.t_frame{j}(1);
             t_e = datarun{2}.block.t_frame{j}(end);
             block.t_sp{j} = t_sp( (t_sp > t_b) & (t_sp < t_e) );
    end


    %%%%%%%%%%%% FOR FULL RASTERS CHANGE HERE !!!! %%%%%%
    if strcmp(stim_vrf,'BW')
        Raster_Blocks = StimulusPars.BW_RasterBlocks;
    elseif strcmp(stim_vrf, 'NSEM')
        Raster_Blocks = StimulusPars.NSEM_RasterBlocks;
    end

    homespikes = zeros( length(Raster_Blocks) , frames*Basepars.spikebins_perstimframe );
    c= 0;
    for k = Raster_Blocks
        c = c+1;
        t_sp_iter = block.t_sp{k} - datarun{2}.block.t_frame{k}(1); % unit of time: sec, 0 for the onset of the block
        for i = 1 : length(t_sp_iter)
            time = t_sp_iter(i);
            bin = ceil ( Basepars.spikebins_perstimframe * time / (Basepars.tstim));
            homespikes ( c , bin) = 1;
        end



    end
    homespikes = homespikes';

    %% MODEL INPUT BUILDING
    %OPT_PARARMS
    p_opt = opt_param.p;
    mu      = p_opt(1);
    mu_bin  = mu * ones ( frames * Basepars.spikebins_perstimframe , 1);
    Stimpars.dt = Basepars.tstim;

    if strcmp(Basepars.k_filtermode, 'rk2')
        kx = filterstimulus_train2AH(p_opt,Basepars,Stimpars,0);
        kx_bin = zeros(frames*Basepars.spikebins_perstimframe , 1);
    end
    if strcmp(Basepars.k_filtermode, 'STA')
        kx = sum(fastconv(Stimpars.movie_ROI,p_opt(Basepars.paramind.L)*Basepars.STA,Basepars.k_spacepixels,frames,Basepars.padval),1)';
        kx = kx';
        kx_bin = zeros(frames*Basepars.spikebins_perstimframe , 1);
    end
    binpf = Basepars.spikebins_perstimframe;
    for i = 1 : length(kx)
        indices = 1 + (i-1)*binpf : i*binpf ;
        kx_bin(indices) = kx(i);
    end
    p_ps     = p_opt(Basepars.paramind.PS);
    psfilter = Basepars.ps_basis * p_ps; 
    ps_gain  = exp(psfilter);

    if Basepars.Coupling
        nNeighbors = length(Basepars.cp_Neighbors);
        neighborspikes = cell( nNeighbors , 1);
        for iNeighbor  = 1 : nNeighbors
            nid = Basepars.cp_Neighbors (iNeighbor) ;
            n_master_idx = find(datarun{1}.cell_ids           == nid);
            n_slave_idx  = find(datarun{2}.cell_ids_map(:,1)  == nid);  
            block.t_sp = cell(StimulusPars.n_blk,1);
            spikes = zeros( length(Raster_Blocks) , frames*Basepars.spikebins_perstimframe );
            if ~isempty(n_slave_idx)
                t_sp     = datarun{2}.spikes{n_slave_idx};        
                for j = 1:StimulusPars.n_blk
                     datarun{2}.block.t_frame{j} = t_frame_interp(datarun{2}.block.trig{j});   %interp stands for interpolation
                     t_b = datarun{2}.block.t_frame{j}(1);
                     t_e = datarun{2}.block.t_frame{j}(end);
                     block.t_sp{j} = t_sp( (t_sp > t_b) & (t_sp < t_e) );
                end
                if strcmp(stim_vrf,'BW')
                    Raster_Blocks = StimulusPars.BW_RasterBlocks;
                elseif strcmp(stim_vrf, 'NSEM')
                    Raster_Blocks = StimulusPars.NSEM_RasterBlocks;
                end


                c = 0 ;
                for k = Raster_Blocks
                    c = c+1;
                    t_sp_iter = block.t_sp{k} - datarun{2}.block.t_frame{k}(1); % unit of time: sec, 0 for the onset of the block
                    for i = 1 : length(t_sp_iter)
                    time = t_sp_iter(i);
                    bin = ceil ( Basepars.spikebins_perstimframe * time / (Basepars.tstim));
                    spikes ( c , bin) = 1;
                    end
                end

                if isfield(Basepars, 'BiDirect_CP') && Basepars.BiDirect_CP
                    shift = Basepars.cp_timebins;
                   spikespart1 = spikes ( :, (shift +1):end);
                   spikespart2 = zeros(length(Raster_Blocks) , shift);
                    spikes      = [spikespart1 , spikespart2 ];
                end
            end
            neighborspikes{iNeighbor}.id     = nid;
            neighborspikes{iNeighbor}.spikes = spikes';

            if isempty(n_slave_idx)
                neighborspikes{iNeighbor}.foundin_vrfset = false;
            end

            %figure; imagesc(spikes);
        end
        microBins_offset = 0;
        CP_Ind  = Basepars.paramind.CP;
        cpno    = Basepars.cp_filternumber;
        for iNeighbor = 1 : nNeighbors
                CP_neigh = CP_Ind ( 1 + (iNeighbor-1) * cpno  : iNeighbor *cpno ); 
                neighborspikes{iNeighbor}.convolved = (grad_basisAH([neighborspikes{iNeighbor}.spikes],Basepars.cp_basis*p_opt(CP_neigh),0));% note: coupling is not examined here.
        end

    end

    %% Simulations

    trials = min(size(homespikes));
    Simulations = cell(trials  , 1);
    simbins = frames * Basepars.spikebins_perstimframe;
    ps_bins = length(ps_gain);
    dt = (1/ Basepars.spikebins_perstimframe) * Stimpars.dt ;
    if ~Basepars.ps_FIX;
        for i_trial = 1 : trials
            lcif_mukx          = mu_bin + kx_bin;
            cpterm = zeros( frames* Basepars.spikebins_perstimframe , 1 );

            if Basepars.Coupling
                for iNeighbor = 1:nNeighbors
                    addin = ( neighborspikes{iNeighbor}.convolved{i_trial} )' ;
                    cpterm = cpterm + addin;
                end
            end
            lcif_mukx_cp = lcif_mukx + cpterm;
            cif_mukxcp0 = exp(lcif_mukx_cp);
            cif_mukxcp_ps = cif_mukxcp0;
            binarycp_simulation = zeros(1,simbins);
            for i = 1 : simbins- ps_bins;
                roll = rand ( 1) ;  
                cif  = cif_mukxcp_ps(i);
                if roll >  exp(-dt*cif);
                    cif_mukxcp_ps(i+1: i + ps_bins) =  cif_mukxcp_ps(i+1: i + ps_bins) .* (ps_gain);
                    binarycp_simulation(i)= 1;

                end
            end
            spikes_cp = size(find(binarycp_simulation));

            Simulations{i_trial}.spikes = binarycp_simulation;
            Simulations{i_trial}.cif    = cif_mukxcp_ps;

            pctdone = i_trial / trials;
            display(sprintf('PercentSimulationsDone_%d',round(100*pctdone)));
        end
    end


    if Basepars.ps_FIX;
            for i_trial = 1 : trials
            lcif_mukx          = mu_bin + kx_bin;
            cpterm = zeros( frames* Basepars.spikebins_perstimframe , 1 );

            if Basepars.Coupling
                for iNeighbor = 1:nNeighbors
                    addin = ( neighborspikes{iNeighbor}.convolved{i_trial} )' ;
                    cpterm = cpterm + addin;
                end
            end

            lcif_mukx_cp = lcif_mukx + cpterm;
            cif_mukxcp0 = exp(lcif_mukx_cp);

            binarycp_simulation = zeros(1,simbins);
            for i = 1 : simbins- ps_bins;
                roll = rand ( 1) ;  
                cif  = cif_mukxcp0(i);
                if roll >  exp(-dt*cif);
                    binarycp_simulation(i)= 1;
                end
            end
            spikes_cp = size(find(binarycp_simulation));
            Simulations{i_trial}.spikes = binarycp_simulation;
            Simulations{i_trial}.cif    = cif_mukxcp0;
            pctdone = i_trial / trials;
            display(sprintf('PercentSimulationsDone_%d',round(100*pctdone)));
            end
    end
%%    PLOTTING MFO
    testseconds     = StimulusPars.nsec_o;
    plotseconds     = 3;
    figs            = floor (testseconds / plotseconds) ;
    plotbins        = floor(Basepars.spikebins_perstimframe *plotseconds /Basepars.tstim );
    for i_fig  = 1 : figs
        figure, clf
        MS = 6;
        bin_start = 1 + (i_fig-1) * plotbins;
        bin_end   = i_fig * plotbins ;   
        sec_start = 0 + (i_fig-1) *plotseconds;
        sec_end   = i_fig *plotseconds;
        x0         = linspace( sec_start , sec_end , plotbins);

        if bin_end > simbins
            bin_end = simbins;
            x0 =  linspace( sec_start , sec_end ,bin_end - bin_start +1);
        end

        %%%%%%%%%%%%
        subplot(211);
        for k = 1:trials
            x = x0 ( find(homespikes(bin_start:bin_end,k)) );
            y = (k/trials) * ones(1, length(x) );
            plot(x,y,'r.','markersize',MS);
           % y = (k/trials)*homespikes(bin_start:bin_end,k);
           % plot(x0,y,'r.','markersize',MS);
            if k == 1, hold on, xlabel('Seconds'); ylabel('RGC Raster'); end
        end
        subplot(212)
        for k = 1:trials
            x = x0 ( find(Simulations{k}.spikes(bin_start:bin_end)) );
            y = (k/trials) * ones(1, length(x) );
            %y = (k/trials)*Simulations{k}.spikes(bin_start:bin_end);
            plot(x,y,'k.','markersize',MS);
            if k == 1, hold on, xlabel('Seconds'); ylabel('GLM'); end
        end
        orient tall
        eval(sprintf('print -dpdf %d_%s_part%d_%s',Basepars.headid,stim_vrf,i_fig,Basepars.fn_save));
    end


end

%{
        %%%% DEFINTIELY WANT TO INHERIT THIS
    max_ind    = ( datarun{1}.vision.sta_fits{master_idx}.mean );
    max_x      = round( max_ind(1)* (movie_width  /StimulusPars.master_width)  );
    max_y      = movie_height - round( max_ind(2)* (movie_height /StimulusPars.master_height) );
    %%%%%%%%%%%%%%%%%%%%%%%
    if rem(GLMPars.K_slen,2) == 1
          offset = (GLMPars.K_slen-1)/2; % off-set from the peak
    else
        offset = GLMPars.K_slen / 2;
    end
    xhigh = (max_x+offset);
    xlow  = (max_x-offset);
    if rem(GLMPars.K_slen,2) == 0
        xlow = xlow +1 ;
    end
    if xhigh > movie_width
           xhigh = xdim; xlow  = xdim - K_slen + 1;
    end
    if xlow < 1
           xlow = 1; xhigh = K_slen ;
    end
    yhigh = (max_y+offset);
    ylow  = (max_y-offset);
    if rem(GLMPars.K_slen,2) == 0
        ylow = ylow +1 ;
    end
    if yhigh > movie_height
           yhigh = ydim; ylow  = ydim - K_slen + 1;
    end
       if ylow < 1
           ylow = 1; yhigh = K_slen ;
       end
    ROI_Xind = (xlow:xhigh);
    ROI_Yind = (ylow:yhigh);
%}
    
    
    
    
    