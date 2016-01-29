%% Bootstrapping script

function results = bootstrap_neuron(basepars,nbootstraps,cellid,nblocks,blocksize,X,save_dir,stim_type,cell_ids,spikes_lng,tstim,gopts)


basepars.maxt = blocksize;
basepars.wholesetsize = basepars.maxt;
basepars.prelimFrames = basepars.maxt;
base_offset = basepars.frame_offset;



blocksdone = [];
blockresults = cell(nblocks,1);


% Figure out the cropping index for this cell
dataidx = find(cell_ids == cellid,1);

basepars.Nneurons = 1;
1;
[stas rstas cpidces box_m box_n ACmat] = crop_sta(spikes_lng{dataidx},X,basepars,struct('dt',tstim),0);
%play_sta(reshape(stas{1},32*64,40),32,64,40);
if (length(cpidces{1}) < basepars.n)
    results = [];
    return;
end
X = X(cpidces{1},:);
basepars.stim_height = sqrt(basepars.n);
basepars.stim_width = sqrt(basepars.n);
basepars.nocropflag = 1;




for k=1:nbootstraps
    
    % Pick a random block
    blocknum = floor(rand*nblocks); % btw 0 and nblocks inclusive
    basepars.frame_offset = base_offset + blocknum*blocksize;
    
    tag = sprintf('cell_%d_bootstrap_%d_block_%dof%d_%s_%s.mat',cellid,k,blocknum+1,nblocks,basepars.filtermode,date);
    
    fprintf('Cell if: %d: Bootstrap run %d: Training on block %d/%d: (frames %d-%d)\n',cellid,k,blocknum+1,nblocks,basepars.frame_offset+1,basepars.frame_offset+basepars.maxt);
    
    oldidx = find(blocksdone == blocknum,1);
    
    if (~isempty(oldidx)) % this block has already been computed - copy the results
        fprintf('This block (%d) has already been computed. Copying old fits.\n',blocknum);
        results.pstar(:,k) = blockresults{oldidx}.pstar;
        results.fstar(k) = blockresults{oldidx}.fstar;
        results.gstar(:,k) = blockresults{oldidx}.gstar;
    else % this is a new block - run the fitting
        blocksdone = [blocksdone;blocknum];
        resultsidx = length(blocksdone);
        [blockresults{resultsidx}.pstar basepars2 blockresults{resultsidx}.fstar blockresults{resultsidx}.gstar] = multiple_neuron(cellid,cellid,basepars,spikes_lng,tstim,cell_ids,stim_type,[],save_dir,gopts,X,tag);
        1;
        results.pstar(:,k) = blockresults{resultsidx}.pstar;
        results.fstar(k) = blockresults{resultsidx}.fstar;
        results.gstar(:,k) = blockresults{resultsidx}.gstar;
        
    end
end

basepars2.cpidces{1} = cpidces{1};
results.basepars = basepars2;
