% NB 2014-10-9
% get data into a reasonable shape to share

for i_exp = 1:3
    for stim_type = {'WN'}
        % Data Type
        GLMType.fit_type = stim_type{1};
        switch i_exp
            case 1
                exp_nm = '2012-08-09-3';
            case 2
                exp_nm = '2012-09-27-3';
            case 3
                exp_nm = '2013-08-19-6';
        end
        
        GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
        GLMType.map_type = 'mapPRJ';
        piece_file = exp_nm(exp_nm ~= '-');
        tstim = .00832750;
        
        %% load basic info
        [StimulusPars, DirPars, datarun_slv, datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
        clear boolean_debug map_type fit_type shead_cellID expname
        inputs.exp_nm    = exp_nm;
        inputs.map_type  = GLMType.map_type;
        inputs.stim_type = GLMType.fit_type;
        inputs.exp_nm       = exp_nm;
        inputs.map_type     = GLMType.map_type;
        DirPars.WN_STAdir   = NSEM_secondaryDirectories('WN_STA', inputs);
        inputs.stim_type    = GLMType.fit_type;
        DirPars.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', inputs);
        clear inputs
        
        %% loading and organizing stimulus
        load(['/Volumes/Lab/Users/Nora/ShareData/Data/CarlosData/' stim_type{1} '-' exp_nm '-StimData.mat'])
        eval(['blockedmoviecell = ' stim_type{1} 'StimData.FitMovie;'])
        eval(['testmovie = ' stim_type{1} 'StimData.TestMovie;'])
        clear WNStimData
        % concatenate the fit movie (different blocks
        n_blocks = length(blockedmoviecell);
        height       = size(blockedmoviecell{1},2);
        width        = size(blockedmoviecell{1},1);
        fitframes    = size(blockedmoviecell{1},3);
        totalframes       = n_blocks * ( fitframes) ;
        concat_fullfitMovie = uint8(zeros(width, height, totalframes)) ;
        for i_blk = 1:n_blocks
            framenums = ( (i_blk -1)*fitframes + 1 ) :  (i_blk *fitframes);
            concat_fullfitMovie(:,:,framenums) = blockedmoviecell{i_blk};
        end
        fitmovie = concat_fullfitMovie;
        clear concat_fullfitMovie blockedmoviecell framenums height i_blk totalframes width
        
        for i_type = 1:2
            if i_type == 1
                files = dir([DirPars.organizedspikesdir '/*OFFPar*']);
                id_start = 8;
            else
                files = dir([DirPars.organizedspikesdir '/*ONPar*']);
                id_start = 7;
            end
            
            %%
            for i_cell = 1:10%length(files)
                
                cell_savename = files(i_cell).name(17:(end-4));
                cid = str2double(cell_savename(id_start:end));
                
                master_idx         = find(datarun_mas.cell_ids == cid);
                stafit_centercoord = ( datarun_mas.vision.sta_fits{master_idx}.mean );
                stafit_sd          = ( datarun_mas.vision.sta_fits{master_idx}.sd   );
                slvdim.height      = StimulusPars.slv.height; slvdim.width = StimulusPars.slv.width;
                [center_coord,~]  = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, StimulusPars.master, slvdim);
                clear master_idx stafit_centercoord slvdim sd stafit_sd stafit_centercoord master_idx
                
                eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, cell_savename));
                %eval(sprintf('load %s/STAandROI_%s.mat STAandROI', DirPars.WN_STAdir, cell_savename));
                spikes = organizedspikes.block.t_sp_withinblock;
                %STA = STAandROI.STA;
                %clear STAandROI
                
                % concatenate spikes
                testspikes = spikes(1:2:end);
                blockedspikes = spikes(2:2:end);
                t_start   = 0;
                T_SP = []; blk_count = 0;
                dur = tstim * fitframes;
                for k = 1:n_blocks
                    blk_count = blk_count + 1;
                    t_sp_full = blockedspikes{k} ; % unit of time: sec, 0 for the onset of the block
                    t_sp      = t_sp_full(t_sp_full >  t_start);
                    t_sp = t_sp - t_start;
                    t_spcontext = t_sp + ( blk_count -1 )*dur;
                    T_SP = [T_SP ; t_spcontext];
                end
                fitspikes = T_SP;
                clear T_SP blk_count blockedspikes k spikes t_sp t_sp_full t_spcontext t_start
                
                %%{
                fittedGLM = glm_fit(fitspikes, fitmovie, [center_coord.x_coord center_coord.y_coord], 'monitor_refresh', 1/tstim);
                fittedGLM.xval = glm_predict(fittedGLM,testmovie, 'testspikes', testspikes);
                temp = corrcoef(conv(sum(fittedGLM.xval.rasters.glm_sim), gausswin(100)),conv(sum(fittedGLM.xval.rasters.recorded), gausswin(100)));
                fittedGLM.xval.corr = temp(2,1);
                clear temp
                close all
                plotfilters(fittedGLM)
                exportfig(gcf, ['/Volumes/Lab/Users/Nora/GLMFits/' piece_file '/' stim_type{1} '/' cell_savename '_filters.eps'], 'Bounds', 'loose', 'Color', 'rgb');
                close all
                plotrasters(fittedGLM.xval, fittedGLM)
                exportfig(gcf, ['/Volumes/Lab/Users/Nora/GLMFits/' piece_file '/' stim_type{1} '/' cell_savename '_rasters.eps'], 'Bounds', 'loose', 'Color', 'rgb');
                close all
                save(['/Volumes/Lab/Users/Nora/GLMFits/' piece_file '/' stim_type{1} '/' cell_savename '.mat'], 'fittedGLM');
                %}
                clear fitspikes center_coord STA cell_savename cid
            end
        end
        clear fitmovie
    end
end

