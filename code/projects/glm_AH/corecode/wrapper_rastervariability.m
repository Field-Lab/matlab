clear; close all; clc
frontoff = 1;
backoff = 1;
for date = 2
    clear datarun_BW; clear datarun_NSEM;
    if date == 1;
        exp_nm = '2012-08-09-3'
        [StimulusPars DirPars datarun_BW]        = Directories_Params_func(exp_nm, 'BW', false);
        [StimulusPars DirPars datarun_NSEM]      = Directories_Params_func(exp_nm, 'NSEM', false);   
    end
    
    if date == 2;
        exp_nm = '2012-09-27-3'
        [StimulusPars DirPars datarun_BW]        = Directories_Params_func(exp_nm, 'BW', false);
        [StimulusPars DirPars datarun_NSEM]      = Directories_Params_func(exp_nm, 'NSEM', false);   

    end
    
	master_Parasols       = union(datarun_BW{1}.cell_types{1}.cell_ids , datarun_BW{1}.cell_types{2}.cell_ids);
    BW_cells_all          = datarun_BW{2}.cell_ids_map(:,1);
    NSEM_cells_all        = datarun_NSEM{2}.cell_ids_map(:,1);
    BWandNSEM             = intersect(BW_cells_all, NSEM_cells_all);
    cells                 = intersect( master_Parasols, BWandNSEM);

    
    for count = 1 : 2
        if count == 1
            stim_vrf = 'BW';
            datarun = datarun_BW;
            nsec_o  = StimulusPars.BW_nsec_o;
        elseif count == 2
            stim_vrf = 'NSEM';
            datarun = datarun_NSEM;
            nsec_o  = StimulusPars.NSEM_nsec_o;
        end
        clear datarun
        [StimulusPars DirPars datarun]      = Directories_Params_func(exp_nm, stim_vrf, false);
        
        
        
        
        for i_cell = cells;
             Cell_ID = i_cell
            isitONParasol = ~isempty( find(datarun{1}.cell_types{1}.cell_ids == Cell_ID) );
            isitOFFParasol = ~isempty( find(datarun{1}.cell_types{2}.cell_ids == Cell_ID) );
            if isitONParasol  && ~isitOFFParasol
                        CTYPE = 'ON-Parasol';
            end
            if ~isitONParasol && isitOFFParasol
                        CTYPE = 'OFF-Parasol';
            end
            output_dir        = DirPars.output_dir; 
            metvar_dirname    = sprintf('%s/%s/%s/metric_variability' ,output_dir , stim_vrf, CTYPE); 
            if (~isdir(metvar_dirname))
                        mkdir(metvar_dirname); 
             end
             if~exist(sprintf('%s/%d_RetRast_PairDist.mat',metvar_dirname,i_cell))
                
           
                plot_logical = true; debug_logical =false;
                [Raster_BinnedLogical, Raster_sptimeSecs ] = load_rasters(Cell_ID, exp_nm, stim_vrf, 10, debug_logical, plot_logical);

              




                trials = max ( size ( Raster_sptimeSecs) );
               %  EVALUATE EACH PAIR METRIC DIST
                   
                RasterCell = Raster_sptimeSecs; raster_time = nsec_o; start_delay = frontoff; end_early = backoff;
                maxpairs =4*trials;
                [Raster_Variability] = raster_variation_metrics( RasterCell, raster_time, start_delay, end_early, maxpairs);
                Raster_Variability.cid       = i_cell;
                Raster_Variability.exp_nm    = exp_nm;
                Raster_Variability.celltype  = CTYPE
                clear RasterCell; clear raster_time; clear start_delay; clear end_early; clear maxpairs;
                display('%%% Computing and Saving Raster Variability Stats %%%')
                save(sprintf('%s/%d_RetRast_PairDist.mat',metvar_dirname,i_cell),'Raster_Variability');
            end
            %end
            %else
            %    load(sprintf('%s/%d_RetRast_PairDist.mat',metvar_dirname,i_cell))
            %    display('%%% Loaded Previously Computed Raster Variability Stats %%%')
            %end
        end
    end
end
            
           