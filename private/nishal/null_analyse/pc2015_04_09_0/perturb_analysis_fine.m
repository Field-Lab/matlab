
addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%% Perturbation ; fine 8, coarse 16 analysis .
% details 
% Fine Stixel 8, Coarse 16
% perturb/base 10,13,14,15,
% null 11
% rgb 7, 12

%% Load movie
interval=2;
rawMovFrames=108000/(2);
[perturbed,height,width,header_size] = get_raw_movie('/Volumes/Lab/Users/bhaishahster/NSbrownian_code_modifiable/rawMovies_final/perturbed_f4c8i1s900.rawMovie',rawMovFrames,1);
perturbed=permute(perturbed,[2,3,1]);
 
[base,height,width,header_size] = get_raw_movie('/Volumes/Lab/Users/bhaishahster/NSbrownian_code_modifiable/rawMovies_final/base_f4c8i1s900.rawMovie',rawMovFrames,1);
base=permute(base,[2,3,1]);

%Repeat frames. 
mov=perturbed;
mov_rep=zeros(size(mov,1),size(mov,2),size(mov,3)*interval);
icnt=0;
for itime=1:size(mov,3)
    if(rem(itime,1000)==1)
    itime
    end
    for iint=1:interval
    icnt=icnt+1;
    
    mov_rep(:,:,icnt)=mov(:,:,itime);
    end
end
perturbed_rep=mov_rep;

%Repeat frames. 
mov=base;
mov_rep=zeros(size(mov,1),size(mov,2),size(mov,3)*interval);
icnt=0;
for itime=1:size(mov,3)
    if(rem(itime,1000)==1)
    itime
    end
    for iint=1:interval
    icnt=icnt+1;
    mov_rep(:,:,icnt)=mov(:,:,itime);
    end
end
base_rep=mov_rep;

clear mov mov_rep


movie = perturbed_rep;
movieLen=size(movie,3);
   movie4D = zeros(size(movie,1),size(movie,2),3,movieLen);
    for iframe=1:size(movie,3);
        if(rem(iframe,1000)==1)
        iframe
        end
        
    movie4D(:,:,1,iframe)=movie(:,:,iframe);
    movie4D(:,:,2,iframe)=movie(:,:,iframe);
    movie4D(:,:,3,iframe)=movie(:,:,iframe);
    end
perturbed_rep_4D = movie4D;

movie = base_rep;
movieLen=size(movie,3);
   movie4D = zeros(size(movie,1),size(movie,2),3,movieLen);
    for iframe=1:size(movie,3);
        if(rem(iframe,1000)==1)
        iframe
        end
        
    movie4D(:,:,1,iframe)=movie(:,:,iframe);
    movie4D(:,:,2,iframe)=movie(:,:,iframe);
    movie4D(:,:,3,iframe)=movie(:,:,iframe);
    end
base_rep_4D = movie4D;

%save('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/movies.mat','base_rep_4D','perturbed_rep_4D','-v7.3');

%%
try
%% Load spikes and format it.
% data003 - perturbed - coarse
mov = (perturbed_rep_4D-127.5)/255;

WN_datafile = '2015-04-09-0/stix8-norefit/data000-from-data000_data003_data004_data005_data008_data009_data016_data017_data006/data000-from-data000_data003_data004_data005_data008_data009_data016_data017_data006';
Null_datafile = '/Volumes/Analysis/2015-04-09-0/stix8-norefit/data003-from-data000_data003_data004_data005_data008_data009_data016_data017_data006';
neuronPath = [Null_datafile,sprintf('/data003-from-data000_data003_data004_data005_data008_data009_data016_data017_data006.neurons')];

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

for ref_cell=1:10%length(InterestingCell_vis_id)
cellID=InterestingCell_vis_id(ref_cell)
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
CellSpkTimes=neuronFile.getSpikeTimes(cellID);
TTL=double(neuronFile.getTTLTimes());

% adjust spikes
spikes=double(CellSpkTimes)/20000;
triggers=TTL/20000;
        spikes_adj=spikes;
        n_block=0;
        for i=1:(length(triggers)-1)
            actual_t_start=triggers(i);
            supposed_t_start=n_block*100/120;
            idx1=spikes > actual_t_start;
            idx2=spikes < triggers(i+1);
            spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
            n_block=n_block+1;
        end
       
% so spikes_adj is the spike times in seconds , adjusted to frame times;
t_bin        = 1/120;
home_sptimes = spikes_adj';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins = home_spbins(find(home_spbins < size(mov,4)) );

   
% home_spbins and mov to calculate STAs.
WNSTA=calc_STA(mov,home_spbins);
    
save(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data003/%d.mat',cellID),'WNSTA','home_spbins');

gmail('bhaishahster@gmail.com', 'A cell done');
end
  gmail('bhaishahster@gmail.com', 'data003 done')

%% Load spikes and format it.
% data005 -base - fine
mov = (base_rep_4D-127.5)/255;

WN_datafile = '2015-04-09-0/stix8-norefit/data000-from-data000_data003_data004_data005_data008_data009_data016_data017_data006/data000-from-data000_data003_data004_data005_data008_data009_data016_data017_data006';
Null_datafile = '/Volumes/Analysis/2015-04-09-0/stix8-norefit/data005-from-data000_data003_data004_data005_data008_data009_data016_data017_data006';
neuronPath = [Null_datafile,sprintf('/data005-from-data000_data003_data004_data005_data008_data009_data016_data017_data006.neurons')];

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

for ref_cell=1:10%length(InterestingCell_vis_id)
cellID=InterestingCell_vis_id(ref_cell)
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
CellSpkTimes=neuronFile.getSpikeTimes(cellID);
TTL=double(neuronFile.getTTLTimes());

% adjust spikes
spikes=double(CellSpkTimes)/20000;
triggers=TTL/20000;
        spikes_adj=spikes;
        n_block=0;
        for i=1:(length(triggers)-1)
            actual_t_start=triggers(i);
            supposed_t_start=n_block*100/120;
            idx1=spikes > actual_t_start;
            idx2=spikes < triggers(i+1);
            spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
            n_block=n_block+1;
        end
       
% so spikes_adj is the spike times in seconds , adjusted to frame times;
t_bin        = 1/120;
home_sptimes = spikes_adj';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins = home_spbins(find(home_spbins < size(mov,4)) );

   
% home_spbins and mov to calculate STAs.
WNSTA=calc_STA(mov,home_spbins);
    
save(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data005/%d.mat',cellID),'WNSTA','home_spbins');


gmail('bhaishahster@gmail.com', 'A cell done');
end
  gmail('bhaishahster@gmail.com', 'data005 done')
 
%% Load spikes and format it.
% data016 - base - fine
mov = (base_rep_4D-127.5)/255;

WN_datafile = '2015-04-09-0/stix8-norefit/data000-from-data000_data003_data004_data005_data008_data009_data016_data017_data006/data000-from-data000_data003_data004_data005_data008_data009_data016_data017_data006';
Null_datafile = '/Volumes/Analysis/2015-04-09-0/stix8-norefit/data016-from-data000_data003_data004_data005_data008_data009_data016_data017_data006';
neuronPath = [Null_datafile,sprintf('/data016-from-data000_data003_data004_data005_data008_data009_data016_data017_data006.neurons')];

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

for ref_cell=1:10%length(InterestingCell_vis_id)
cellID=InterestingCell_vis_id(ref_cell)
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
CellSpkTimes=neuronFile.getSpikeTimes(cellID);
TTL=double(neuronFile.getTTLTimes());

% adjust spikes
spikes=double(CellSpkTimes)/20000;
triggers=TTL/20000;
        spikes_adj=spikes;
        n_block=0;
        for i=1:(length(triggers)-1)
            actual_t_start=triggers(i);
            supposed_t_start=n_block*100/120;
            idx1=spikes > actual_t_start;
            idx2=spikes < triggers(i+1);
            spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
            n_block=n_block+1;
        end
       
% so spikes_adj is the spike times in seconds , adjusted to frame times;
t_bin        = 1/120;
home_sptimes = spikes_adj';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins = home_spbins(find(home_spbins < size(mov,4)) );

   
% home_spbins and mov to calculate STAs.
WNSTA=calc_STA(mov,home_spbins);
    
save(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data016/%d.mat',cellID),'WNSTA','home_spbins');


gmail('bhaishahster@gmail.com', 'A cell done');
end
  gmail('bhaishahster@gmail.com', 'data016 done')


%% Load spikes and format it.
% data017 - base - coarse
mov = (perturbed_rep_4D-127.5)/255;

WN_datafile = '2015-04-09-0/stix8-norefit/data000-from-data000_data003_data004_data005_data008_data009_data016_data017_data006/data000-from-data000_data003_data004_data005_data008_data009_data016_data017_data006';
Null_datafile = '/Volumes/Analysis/2015-04-09-0/stix8-norefit/data017-from-data000_data003_data004_data005_data008_data009_data016_data017_data006';
neuronPath = [Null_datafile,sprintf('/data017-from-data000_data003_data004_data005_data008_data009_data016_data017_data006.neurons')];

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end

for ref_cell=1:10%length(InterestingCell_vis_id)
cellID=InterestingCell_vis_id(ref_cell)
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
CellSpkTimes=neuronFile.getSpikeTimes(cellID);
TTL=double(neuronFile.getTTLTimes());

% adjust spikes
spikes=double(CellSpkTimes)/20000;
triggers=TTL/20000;
        spikes_adj=spikes;
        n_block=0;
        for i=1:(length(triggers)-1)
            actual_t_start=triggers(i);
            supposed_t_start=n_block*100/120;
            idx1=spikes > actual_t_start;
            idx2=spikes < triggers(i+1);
            spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
            n_block=n_block+1;
        end
       
% so spikes_adj is the spike times in seconds , adjusted to frame times;
t_bin        = 1/120;
home_sptimes = spikes_adj';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins = home_spbins(find(home_spbins < size(mov,4)) );

   
% home_spbins and mov to calculate STAs.
WNSTA=calc_STA(mov,home_spbins);
    
save(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data017/%d.mat',cellID),'WNSTA','home_spbins');

gmail('bhaishahster@gmail.com', 'A cell done');

end
  gmail('bhaishahster@gmail.com', 'data017 done')

catch
     gmail('bhaishahster@gmail.com', 'Error in code')
end

%%
cellID = 1022;
pert_offset0 = load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data003/%d.mat',cellID));

pert_offset2 = load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data005/%d.mat',cellID));

base_offset0 = load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data016/%d.mat',cellID));

base_offset2 = load(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_04_09_0/data017/%d.mat',cellID));

STAmagnitude = sum(base_offset0.WNSTA.^2,3);
max_val = max(STAmagnitude(:));
[r,c]=find(max_val==STAmagnitude);
mask = 0*STAmagnitude;
mask(r(1)-20:r(1)+20,c(1)-20:c(1)+20)=1;
mask=logical(repmat(mask,[1,1,30]));

mask_off = 0*STAmagnitude;
mask_off(r(1)-20:r(1)+20,c(1)-20:c(1)+20)=1; % Doubt
mask_off=logical(repmat(mask_off,[1,1,30]));


[v,STAidx]=max(abs(squeeze(base_offset0.WNSTA(r(1),c(1),:))));

figure;
subplot(3,2,1);
a=pert_offset0;
b=reshape(a.WNSTA(mask),[41,41,30]);
imagesc(b(:,:,STAidx));
caxis([min(a.WNSTA(:)),max(a.WNSTA(:))]);
colormap gray
colorbar

subplot(3,2,2);
a=pert_offset2;
b=reshape(a.WNSTA(mask_off),[41,41,30]);
imagesc(b(:,:,STAidx));
caxis([min(a.WNSTA(:)),max(a.WNSTA(:))]);
colormap gray
colorbar

subplot(3,2,3);
a=base_offset0;
b=reshape(a.WNSTA(mask),[41,41,30]);
imagesc(b(:,:,STAidx));
caxis([min(a.WNSTA(:)),max(a.WNSTA(:))]);
colormap gray
colorbar

subplot(3,2,4);
a=base_offset2;
b=reshape(a.WNSTA(mask_off),[41,41,30]);
imagesc(b(:,:,STAidx));
caxis([min(a.WNSTA(:)),max(a.WNSTA(:))]);
colormap gray
colorbar


subplot(3,2,5);
xx =  base_offset0.WNSTA(mask)/norm( base_offset0.WNSTA(mask));
yy = pert_offset0.WNSTA(mask);
STAdiff=yy - (yy'*xx)*xx;  % adjust scales
b=reshape(STAdiff,[41,41,30]);
imagesc(b(:,:,STAidx));
caxis([min(STA(:)),max(STA(:))]);
colormap gray
colorbar

subplot(3,2,6);
xx =  base_offset2.WNSTA(mask)/norm( base_offset2.WNSTA(mask));
yy = pert_offset2.WNSTA(mask);
STAdiff=yy - (yy'*xx)*xx;  % adjust scales
b=reshape(STAdiff,[41,41,30]);
imagesc(b(:,:,STAidx));
caxis([min(STA(:)),max(STA(:))]);
colormap gray
colorbar
