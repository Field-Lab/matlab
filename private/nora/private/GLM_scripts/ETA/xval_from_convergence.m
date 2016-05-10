function [movie, xval, res] = xval_from_convergence
Conv = 20;
Cutoff = -1; % BPS cutoff
datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_noCP_p8IDp8/standardparams/NSEM_mapPRJ/'; %Mod1Max1e4
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
n_cells = 15;
res.cells = zeros(n_cells, 1);
res.centers = zeros(n_cells, 2);
i_cell = 0;

for exp = 1
    matfiles = dir([datapath exp_names(exp,:) 'conv_blocks_' num2str(Conv) '/*.mat']);
    load([datapath exp_names(exp,:) 'conv_blocks_' num2str(Conv) '/' matfiles(1).name]);
        
    % Load and process stimulus
    exp_nm = exp_names(exp,1:(end-1));
    [StimulusPars, exp_info] = StimulusParams(exp_nm, 'NSEM', fittedGLM.GLMType.map_type);
    [blockedmoviecell, inputstats, origmatfile] = loadmoviematfile(exp_nm , 'NSEM', fittedGLM.GLMType.cone_model,'fitmovie');
    % StimulusPars.slv.FitBlocks = StimulusPars.slv.FitBlocks((Conv+1):57);
    movie = blockedmoviecell(StimulusPars.slv.FitBlocks((Conv+1):57)/2);
    
    % Find where the spikes are
    secondDir.exp_nm    = exp_nm;
    secondDir.map_type  = fittedGLM.GLMType.map_type;
    secondDir.stim_type = 'WN';
    secondDir.fitname   = fittedGLM.GLMType.fitname;
    Dirs.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir);
    
    for file = 1%:length(matfiles)
        load([datapath exp_names(exp,:) 'conv_blocks_' num2str(Conv) '/' matfiles(file).name]);
        
        if fittedGLM.xvalperformance.glm_normedbits > Cutoff  % Only use the fit if it has relatively converged
            i_cell = i_cell + 1
            res.cells(i_cell) = fittedGLM.cellinfo.cid;
            res.centers(i_cell,1) = fittedGLM.cellinfo.slave_centercoord.y_coord;
            res.centers(i_cell,2) = fittedGLM.cellinfo.slave_centercoord.x_coord;
            
            if fittedGLM.GLMType.CouplingFilters
                for i = 1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
                    eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir, fittedGLM.cellinfo.pair_savename{i})); 
                    neighbor_spikes(i,:) = organizedspikes.block.t_sp_withinblock(StimulusPars.slv.FitBlocks((Conv+1):57));
                end
            else 
                neighbor_spikes = 0;
            end
            
            eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir, fittedGLM.cell_savename));
            spikes = organizedspikes.block.t_sp_withinblock(StimulusPars.slv.FitBlocks((Conv+1):57));
            fittedGLM.center_coord = fittedGLM.cellinfo.slave_centercoord;
            fittedGLM.inputstats.mu_avgIperpix = 64;
            fittedGLM.inputstats.range = 255;
            for block = 1:length(movie)
                xval{i_cell,block} = glm_predict(fittedGLM, movie{block}.matrix, 'testspikes', spikes(block), 'trials', 1, 'predict', false);%, 'neighborspikes', neighbor_spikes(:,block));
            end
        end
        
    end
    
end
end

function concat_fitmovie   = subR_concat_fitmovie_fromblockedcell(blockedmoviecell , FitPars)
% AKHeitman 2014-04-14
% Concatenate the fit movie (different blocks)
% FitPars needs
%   .width
%   .height
%   .FitBlocks
%   .novelblocks
%   .fitframes

height       = FitPars.height;
width        = FitPars.width;
fitblocks    = FitPars.FitBlocks;
fitframes    = FitPars.fitframes;
novelblocks  = FitPars.NovelBlocks;

fitframesperblock = length(fitframes) ;
totalframes       = length(fitblocks) * ( fitframesperblock) ;
concat_fullfitMovie = uint8(zeros(width, height, totalframes)) ;
for i_blk = fitblocks
        blkind = find(fitblocks == i_blk);
        framenums = ( (blkind -1)*fitframesperblock + 1 ) :  (blkind *fitframesperblock);  
        n_blkind = find(novelblocks == i_blk);
        concat_fullfitMovie(:,:,framenums) = blockedmoviecell{n_blkind}.matrix (:,:, fitframes);    
end

concat_fitmovie = concat_fullfitMovie;

end