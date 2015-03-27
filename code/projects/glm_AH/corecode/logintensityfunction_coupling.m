function [lcif_cp, lcif_perneighbor] = logintensityfunction_coupling (pstar, Basepars, rastertype_string, boolean_debug)

% Checked 2012-12-07   AK Heitman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs: pstar (final param values of the GLM)
%        datarun  in the standard cell form {1} {2} master slve (spike
%              timing)
%        Coulping Basis
%        Raster_blocks indicating which blocks to use to build the raster        
%        spikebins    net number of bins to be used
%        block_duratoin   time in seconds of the raster
%        logical_BiDirectCP    whether or not the coupling terms were
%                   optional argument .. default is false!!
%Calls:  grad_basisAH
%Outputs: lcif_cp  lcif_perneighbor
%  lcif_cp:  Cell with cell num equal to raster trials
%            {trial}.couplingterm  scalar function of spike bins
%                   all inclusive coupling term fed into GLM            
%            {trial}.blocknum     the block number in terms of datarun
%  lcif_perneighbor:
%            large structure with spike times and individual filters etc.


% Flow of Operations
% 1. create binned spike structure for all the neighbors
% 2. convolved spikes with corresponding filter  Calling !!grad_basisAH!!
% 3. sum all neighbors to get scalar function of spike bins    lcif_cp



% Load up datarun .. set some variables
[Slv_StimPars DirPars datarun]  = Directories_Params_func_all(Basepars.exp_nm, boolean_debug, rastertype_string);
spikebins_perblock = Basepars.spikebins_perstimframe * (Slv_StimPars.ntb_o-1) * Slv_StimPars.frames_pertrigger;
block_duration     = Slv_StimPars.nsec_o;
Paramind  = Basepars.paramind; Neighbor_IDs = Basepars.cp_Neighbors  ; %% should be all there.. all misses should be eliminated
Coupling_Basis = Basepars.cp_basis;
Raster_Blocks = Slv_StimPars.RasterBlocks;        
% Setting Parameters of the program
cp_params      = pstar(Paramind.CP);
totalblocks    = max(size(datarun{2}.block.trig));
nNeighbors     = length(Neighbor_IDs);
neighborspikes = cell( nNeighbors , 1);

% Make sure we use correct type of Coupling
%if ~exist('logical_BiDirectCP', 'var')
 %   logical_BiDirectCP = false;
%end


if ~isfield(Basepars, 'BiDirect_CP')
    logical_BiDirectCP = false;
else
    logical_BiDirectCP = Basepars.BiDirect_CP;
end

% Binning the Neighbor Spikes from the datarun
for i_Neighbor  = 1 : nNeighbors
	nid = Neighbor_IDs(i_Neighbor) ;
	n_master_idx = find(datarun{1}.cell_ids           == nid);
	n_slave_idx  = find(datarun{2}.cell_ids_map(:,1)  == nid);  
    
    
    isitONParasol  = ~isempty( find(datarun{1}.cell_types{1}.cell_ids == nid) );
    isitOFFParasol = ~isempty( find(datarun{1}.cell_types{2}.cell_ids == nid) );
    if isitONParasol  && ~isitOFFParasol
        CTYPE = 'ON-Parasol';
    end
    if ~isitONParasol && isitOFFParasol
        CTYPE = 'OFF-Parasol';
    end
    
    spikes       = false( length(Raster_Blocks) , spikebins_perblock);
    sptimes_secs = cell(length(Raster_Blocks), 1);
	if ~isempty(n_slave_idx)
        t_sp     = datarun{2}.spikes{n_slave_idx};  
        block.t_sp = cell(totalblocks,1);
        spikes       = false( length(Raster_Blocks) , spikebins_perblock);
        sptimes_secs = cell(length(Raster_Blocks), 1);
        
        
        for j = 1:totalblocks
                     datarun{2}.block.t_frame{j} = t_frame_interp(datarun{2}.block.trig{j});   %interp stands for interpolation
                     t_b = datarun{2}.block.t_frame{j}(1);
                     t_e = datarun{2}.block.t_frame{j}(end);
                     block.t_sp{j} = t_sp( (t_sp > t_b) & (t_sp < t_e) );
        end
        c = 0 ;
        for k = Raster_Blocks
            c = c+1;
            t_sp_iter = block.t_sp{k} - datarun{2}.block.t_frame{k}(1); % unit of time: sec, 0 for the onset of the block
            for i = 1 : length(t_sp_iter)
                    time = t_sp_iter(i);
                    bin = ceil( spikebins_perblock *  (time / block_duration) );
                    spikes ( c , bin) = true;
                    
            end
            sptimes_secs{c}.times = t_sp_iter;
            sptimes_secs{c}.block = k;
        end

        if logical_BiDirectCP 
            display('%%%%%%%%%%%%%%%%adjustmentinSimforBiDirectCoupling%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') 
            shift       = floor(size(Coupling_Basis,1)/2);
            spikespart1 = spikes ( :, (shift +1):end);
            spikespart2 = zeros(length(Raster_Blocks) , shift);
            spikes      = [spikespart1 , spikespart2 ];
        end
    end
	neighborspikes{i_Neighbor}.id               =  nid;
    neighborspikes{i_Neighbor}.celltype         = CTYPE;
	neighborspikes{i_Neighbor}.binned_spikes    = spikes';
    neighborspikes{i_Neighbor}.bins             = spikebins_perblock;
    neighborspikes{i_Neighbor}.raster_blocknum  = Raster_Blocks;
    neighborspikes{i_Neighbor}.spiketimes_secs  = sptimes_secs;
    neighborspikes{i_Neighbor}.foundin_vrfset   = true;
    

	if isempty(n_slave_idx)
        neighborspikes{i_Neighbor}.foundin_vrfset = false;
	end

            %figure; imagesc(spikes);
end

% Conolving the logicaly binned spikes with the cp_params*BasidFunction
microBins_offset = 0;
cpno    = size(Coupling_Basis,2);
for i_Neighbor = 1 : nNeighbors
                singleneighbor = (i_Neighbor-1)*cpno + 1 : i_Neighbor*cpno;
                cp_parsingle   = cp_params(singleneighbor);
                cp_filter      = Coupling_Basis*cp_parsingle;
                neighborspikes{i_Neighbor}.convolved_lcifterm = (grad_basisAH([neighborspikes{i_Neighbor}.binned_spikes],cp_filter,0))';% note: coupling is not examined here.
                neighborspikes{i_Neighbor}.cpfilter           =  cp_filter;
end
        

% Summing to give final scalar output as function of spike bins time
% lcif_cp  the all inclusive coupling term!
clear lcif_cp
lcif_cp  = cell( length(Raster_Blocks), 1); 
for i_trial = 1 : length(Raster_Blocks)
    cpterm = zeros(spikebins_perblock , 1);
    
	for i_Neighbor = 1:nNeighbors
        addin = ( neighborspikes{i_Neighbor}.convolved_lcifterm{i_trial}') ;
        cpterm = cpterm + addin;
    end
    lcif_cp{i_trial}.couplingterm    = cpterm;
    lcif_cp{i_trial}.blocknum = Raster_Blocks(i_trial);
   
end

lcif_perneighbor = neighborspikes;