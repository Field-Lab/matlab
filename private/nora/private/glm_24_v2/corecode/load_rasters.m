function [Raster_BinnedLogical, Raster_sptimeSecs ] = load_rasters(Cell_ID, expnm_string, rastertype_string, bins_perstimframe, logical_debug, plot_logical);
% Complte Checked% AK Heitman   11-29-2012

%INPUTS:
%rastertype_string  =   either 'BW' or 'NSEM'
%expnm_string         something of fomr '2012-09-27-3'
%bins_perstimframe   important parameter.. should be small enough so that
  %at most one spike occurs per bin (>= 5 atleast)
   % in GLM context equals  Basepars.spikebins_perstimframe


%OPERATIONS DONE:   
%Loads datarun
%gives block structure to inidivudal spikes 
%uses Slv_StimPars.RasterBlocks associated from Directories_Params_func_all

% Will load up a Logical Matrix .. each row is a trial each column a bin
% Will also load a cell whose element are each the sptimes in secs for different trials


%%%%%%%%%%%%%%%%%%%%%%%%%%%   LOAD DATARUN AND SOME PARAMS
cid = Cell_ID;
if ~exist('logical_debug', 'var')
    logical_debug = false;
end
if ~exist('plot_logical','var')
    plot_logical = true;
end
[Slv_StimPars , DirPars, datarun] = Directories_Params_func_all(expnm_string, logical_debug, rastertype_string);  
totalblocks                           = Slv_StimPars.n_blk;
t_stim                                = Slv_StimPars.tstim;
%if strcmp(rastertype_string ,'BW')
%    rasterblocks = Slv_StimPars.BW_RasterBlocks; 
%    nsec_o = Slv_StimPars.BW_nsec_o; nsec_e = Slv_StimPars.BW_nsec_e; 
%    ntb_o  = Slv_StimPars.BW_ntb_o;  ntb_e = Slv_StimPars.BW_ntb_e; 
%    ntb_oe = ntb_o + ntb_e; nsec_oe = nsec_o + nsec_e;   
%end
%if strcmp(rastertype_string ,'NSEM')
%    rasterblocks = Slv_StimPars.NSEM_RasterBlocks; 
%    nsec_o = Slv_StimPars.NSEM_nsec_o; nsec_e = Slv_StimPars.NSEM_nsec_e; 
%    ntb_o  = Slv_StimPars.NSEM_ntb_o;  ntb_e = Slv_StimPars.NSEM_ntb_e; 
%    ntb_oe = ntb_o + ntb_e; nsec_oe = nsec_o + nsec_e;    
%end

slave_idx    =  find(datarun{2}.cell_ids_map(:,1) == cid);
master_idx   =  find(datarun{1}.cell_ids          == cid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%   GIVE BLOCK STRUCTURE TO SPIKES
block.t_sp = cell(totalblocks,1);
t_sp     = datarun{2}.spikes{slave_idx};
for j = 1:totalblocks
	datarun{2}.block.t_frame{j} = t_frame_interp(datarun{2}.block.trig{j});   %interp stands for interpolation
	t_b = datarun{2}.block.t_frame{j}(1);
	t_e = datarun{2}.block.t_frame{j}(end);
	block.t_sp{j} = t_sp( (t_sp > t_b) & (t_sp < t_e) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%   CREATE BOTH RASTERS (SPTIMES AND LOGICAL BINS) 
sptimeraster = cell(length(Slv_StimPars.RasterBlocks),1);
bin_num      = (Slv_StimPars.ntb_o-1) * Slv_StimPars.frames_pertrigger * bins_perstimframe;
homespikes   = false( length(Slv_StimPars.RasterBlocks) , ceil(bin_num) );
c= 0;
for k = Slv_StimPars.RasterBlocks
        c = c+1;
        t_sp_iter = block.t_sp{k} - datarun{2}.block.t_frame{k}(1); % unit of time: sec, 0 for the onset of the block
        sptimeraster{c}.block      =   Slv_StimPars.RasterBlocks(c);
        sptimeraster{c}.spikes     = t_sp_iter; 
        sptimeraster{c}.spikecount = length(t_sp_iter);
        sptimeraster{c}.blocklength= Slv_StimPars.nsec_o; 
        sptimeraster{c}.sprate     = length(t_sp_iter) / Slv_StimPars.nsec_o;
        for i = 1 : length(t_sp_iter)
            time = t_sp_iter(i);
            bin = ceil((time/Slv_StimPars.nsec_o) * bin_num);
            homespikes ( c , bin) = true;
        end
        sptimeraster{c}.binned_logical = homespikes(c,:);
end
homespikes = homespikes';


%%%%%%%%%%%%%%%%%%%%%%%%%%%   DIAGNOSTIC PLOT TO SEE RASTER IS WORKING
if plot_logical
    figure
    max_sec = min(5, Slv_StimPars.nsec_o);
    MS = 6; hold on;  xlabel('Seconds'); ylabel('Trial Number'); title('Diagnostic Raster')
    for i = 1 : length(Slv_StimPars.RasterBlocks);
        sptimes = sptimeraster{i}.spikes (find ( sptimeraster{i}.spikes < max_sec) );
        y_base  = ones(size(sptimes));
        plot(sptimes, i * y_base , 'r.', 'markersize',MS);
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%   ASSIGN TO OUTPUT
Raster_BinnedLogical = homespikes';
Raster_sptimeSecs    = sptimeraster;

