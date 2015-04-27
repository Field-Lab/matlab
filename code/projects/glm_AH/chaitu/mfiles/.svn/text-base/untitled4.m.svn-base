% Look at pairwise log-likelihoods

ncells = 6;
nconditions = 4;

datasetnames = fliplr({'X04', 'X05', 'X06', 'X07'});
spikelistnames = fliplr({'spikes09', 'spikes10', 'spikes11', 'spikes13'});
basepath = '/lcv/data/chaitu/glm/single_neuron/2001_dataset/';

ntestframes = 1000;


LGP = zeros(nconditions,nconditions,ncells);
LGPmix = zeros(nconditions,ncells);

for i=1:ncells

    for j=1:nconditions

        datapath = strcat(basepath,'wn0',num2str(nconditions-j+1+3),'_fits/');
        fitlist = what(datapath); fitlist = fitlist.mat;        
        idx = strmatch(strcat('cells_',num2str(i)),fitlist);
        datafile = strcat(datapath,fitlist{idx})
        
        load(datafile);
        
        for k=1:nconditions
        
            X = evalin('base',datasetnames{k});
            spcell = cell(1,1);
            sp = evalin('base',spikelistnames{k});
            spcell{1} = sp{i};
            basepars.crop_idx = [1];
            [LGP(j,k,i) cifs kx lg_p_sep] = eval_glm_ll(master_pstar,basepars,X(:,1:ntestframes),spcell,stimpars.dt,0);
            
            
        end
        
        datapath = strcat(basepath,'wnmix_fits/50K_frames/');
        fitlist = what(datapath); fitlist = fitlist.mat;
        idx = strmatch(strcat('cells_',num2str(i)),fitlist);
        datafile = strcat(datapath,fitlist{idx})
        
        load(datafile);
        X = evalin('base',datasetnames{j});
        spcell = cell(1,1);
        sp = evalin('base',spikelistnames{j});
        spcell{1} = sp{i};
        basepars.crop_idx = [1];
        [LGPmix(j,i) cifs kx lg_p_sep] = eval_glm_ll(master_pstar,basepars,X(:,1:ntestframes),spcell,stimpars.dt,0);
        
    end
        
        
    
    
    
end
