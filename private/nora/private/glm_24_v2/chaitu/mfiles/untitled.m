sep_dir = '/lcv/data/chaitu/glm/single_neuron/single_05_22_2009/1f/';
nonsep_dir = '/lcv/data/chaitu/glm/single_neuron/06_22_2009_nonsep/1f/';

Asep = what(sep_dir);
Asep = Asep.mat;

Anonsep = what(nonsep_dir);
Anonsep = Anonsep.mat;

nsets = min(length(Anonsep),length(Asep));


Fsep = [];
Fnonsep = [];
cellids = [];

% Load stimulus and spike data
load /lcv/data/chaitu/1f_72Kframes_08_2008.mat
load ~/data/shlens_newdata/2008-08-27-6-Spikes.mat;

for j=1:length(Anonsep)

    str_nonsep = Anonsep{j}
    
    dividers = find(str_nonsep == '_');
    
    cellid = str2num(str_nonsep(dividers(1)+1:dividers(2)-1));
    
    sepidx = strmatch(str_nonsep(1:dividers(2)),Asep)
    
    if (isempty(sepidx))
        continue;
    end
    
    str_sep = Asep{sepidx(1)};
    load(strcat(sep_dir,str_sep));

    Fsep = [Fsep; -fstar/basepars.maxt];
    [lg_p_tst D modelD modelPSTH PSTH uisis_tst urisis_tst kxvar_tst] = evaluate_model(X,pstar,basepars,stimpars,trainpars,spikes_lng,spikes_rpt,...
                                                                        psth_dt,nrepeatframes,fulltrial_duration,nframes_correction,...
                                                                        ntrials,ntrials_glm,cellids,cell_ids)    
    
    load(strcat(nonsep_dir,str_nonsep));
    Fnonsep = [Fnonsep; -fstar/basepars.maxt];        
    
    cellids = [cellids; cellid];
    
end
    
%%

ncells = length(cellids);

for j=1:ncells
    
    % Separable model
    load(strcat(sep_dir,str_sep
    
    