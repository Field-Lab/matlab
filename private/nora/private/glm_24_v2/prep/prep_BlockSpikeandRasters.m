% function  cd    prep_BlockSpikeandRasters  
%               2013-12-04  AKHeitman (previously EDOI) 
%
%               
%
% usage:        prep_BlockSpikeandRasters(string_date, slv_type, map_type)
%               plot Rasters and load Spike Times on per Cell basis
%               saves blocked spikes per cell into structure
%               "organizedspikes"
%
% arguments:    slv_type is either 'BW' or 'NSEM'
%               map_type is either 'mapEI' or 'mapPRJ'
%
%
% arguments:    Directories_Params_v18_split
%               func_ras_blkdsnAH
%
% outputs:      NSEM_Projects/GLM/2012-08-09-3/BW_mapEI/BlockSpikesandRasters
%               organizedspikes_celltype_cellid%.mat               
%
% paths:        run glmpath_18.m before everything
%
%
% Previous:     prep_RetinaRasters (glm_AH_18a, glm_AH_17)
% EDOI, 2012-02-09.
% AKHEITMAN 2012-09-11 Clean Up

%%%%% READ %%%%%
% First pass run through the datarun collection.  Verify that we properly
% accounted for triggers and block structure by plotting the data.  Need to
% make sure that the rasters appear correctly. Organized spikes, just puts all the spikes into their  

% Organizedspikes.mat orgnanizes spikes by block
%  - also a general repository for information about the cell

%{
% silly path stuffdataru
% path('/netapp/snle/lab/temp_matlabcodeAH/glm_AH_18', path)
% path('/netapp/snle/lab/temp_matlabcodeAH/glm_AH_18/prep', path)
% path('/netapp/snle/lab/temp_matlabcodeAH/glm_AH_18/bookkeep',path)
% path('/netapp/snle/lab/temp_matlabcodeAH/glm_AH_18/corecode',path)
% path('/netapp/snle/lab/temp_matlabcodeAH/generalcomputations_AKH',path)
%}
%%%%% Notes %%%%%%%%%
% USE MASTER CLASSIFICATION TO CLASSIFY SLAVE DATA 
% SLAVE DATA IS WHERE WE FIT THE MODEL AND TEST IT ON REPEATS
% LOAD UP THE SPIKE TIMES OF THE ENSLAVED DATA AND PLOT RASTERS
% CELL TYPES ARE SPLIT UP
function prep_BlockSpikeandRasters(string_date, slv_type, map_type)

%exp_nm = '2012-09-27-3'; fit_type = 'BW'; map_type = 'mapPRJ'; 
exp_nm = string_date;
boolean_debug = false;
% CTYPE = {'On Parasol','Off Parasol','On Midget','Off Midget', 'SBC', 'nc5', 'nc4', 'nc6', 'Unclassified'}; %
CTYPE = {'On Large', 'Off Large', 'nc5', 'nc4', 'nc6', 'Unclassified'}; %

% CTYPE_sname = { 'ONPar' , 'OFFPar','ONMid','OFFMid','SBC', 'nc5', 'nc6', 'Unknown'};  % short name
CTYPE_sname = {'ONLarge', 'OFFLarge', 'nc5','nc4', 'nc6', 'Unknown'};  % short name

%%% Just or debugging
%{ 
just for debuuging
exp_nm = '2012-08-09-3'; fit_type = 'BW'; map_type = 'mapPRJ';
shortlist_cellID= {[1276], [1471], [7652], [648]}; % [271];
CTYPE = {'On-Parasol','Off-Parasol','On-Midget','Off-Midget'}; %
CTYPE_sname = { 'ONPar' , 'OFFPar','ONMid','OFFMid' };  % short name
%}
[StimulusPars, DirPars, datarun_slv, datarun_mas] = Directories_Params_v19_split(exp_nm, slv_type, map_type)
if strcmp(slv_type,'BW')
    Stim = StimulusPars.BW;
elseif strcmp(slv_type, 'NSEM')
    Stim = StimulusPars.NSEM;
end

d_save = sprintf('%s/BlockSpikesandRasters',DirPars.output_dir);
if ~exist(d_save,'dir')
       mkdir(d_save);
end
rastdir = sprintf('%s/rasters',d_save);
if ~exist(rastdir,'dir')
       mkdir(rastdir);
end

%%

%%% This core computation should already be portable
%%% SAVE THE INFO IN ONE TESTRUN STRUCTURE PER CELL; PLOT RASTERS
%%% OUTPUT GOES INTO snl-e/glm/output/2012-08-21-1/BW
for n_ctype = 1:length(CTYPE)
   
   [celltype_index, CID] = celltype_id_AH(CTYPE{n_ctype}, datarun_mas{1}.cell_types);
    
    
 %  CID = datarun_mas.cell_types{n_ctype}.cell_ids;
   if exist('shortlist_cellID', 'var')
       CID = shortlist_cellID{n_ctype};
   end
   
   %% CYCLE THROUGH ALL CELLS OF THE TPE
   for n_cid = 1:length(CID)
       try
      %% Load data by cell and save ... sort spikes into blocks
      % Also repository for as much information about the cell as possible
      clear organizedspikes
      organizedspikes.cell_id       = CID(n_cid);
      organizedspikes.cell_type     = CTYPE{n_ctype}
      organizedspikes.savename      = sprintf('%s_%d',CTYPE_sname{n_ctype},CID(n_cid));
      organizedspikes.maptype       = map_type;
      organizedspikes.maptype_note  = 'Celltype defined in a a seperate master run .. transferred either by EI mapping or projections mapping  (mapEI or mapPRJ)';
      organizedspikes.no_trig_bl    = [Stim.ntb_o, Stim.ntb_e];
      organizedspikes.notrigbl_note = 'Number of triggers in the alternating blocks of  raster and fit';
      organizedspikes.t_trig        = datarun_slv.triggers;
      organizedspikes.cell_idx      = find(datarun_slv.cell_ids==organizedspikes.cell_id);
      organizedspikes.cell_idx_note = 'index within datarun.cell_ids of cell';
      organizedspikes.t_sp          = datarun_slv.spikes{organizedspikes.cell_idx};
      organizedspikes.t_sp_note     = 'raw spike times in seconds';
      organizedspikes.block.t_sp_note               = 'raw spike times divided by block';
      organizedspikes.block.t_sp_withinblock_note   = 'spike times within blocks, where beginning of each block is time 0';
      organizedspikes.block.t_sp = cell(Stim.n_blk,1); 
      organizedspikes.block.t_sp_withinblock = cell(Stim.n_blk,1); 
      for i_block = 1:Stim.n_blk
          clear t_b t_e  
         t_b = datarun_slv.block.t_frame{i_block}(1);
         t_e = datarun_slv.block.t_frame{i_block}(end);
         organizedspikes.block.t_sp{i_block} = organizedspikes.t_sp( (organizedspikes.t_sp > t_b) & (organizedspikes.t_sp < t_e) );
         organizedspikes.block.t_sp_withinblock{i_block} = organizedspikes.block.t_sp{i_block} - t_b; 
      end
      organizedspikes.block.t_frame = datarun_slv.block.t_frame;
      
      
      %% SAVE AND PRINT OUT RASTERS  (STATIC, PSTH, NOVEL, FIRING RATES)
      eval( sprintf( 'save %s/organizedspikes_%s organizedspikes' , d_save , organizedspikes.savename ) );
      pdf_name = sprintf('%s/rasters/rast_%s',d_save,organizedspikes.savename);
      
      SPars = Stim;
      hack_printRasters(organizedspikes,SPars,pdf_name)
       catch
           disp('missing cell id')
       end
       
   end

end


end



