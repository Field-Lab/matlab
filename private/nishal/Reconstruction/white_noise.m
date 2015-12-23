%% For Linear reconstruction

% Load data and STAs
startup_null_analyse_tenessee
%% Load data


datafile = '2012-08-09-3/data002'%'2012-08-09-3/data002'; %'2012-08-21-1/data003';

opt=struct('load_all',true);
datarun=load_data(datafile,opt)
datarun= load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun=load_params(datarun)
datarun = set_polarities(datarun);

%get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
%STAs!
cellType=datarun.cell_types{1}.name;% 'all'
matlab_cell_ids=get_cell_indices(datarun,cellType);
%stas=datarun.stas.stas(matlab_cell_ids);
%n_cell=length(stas);



%% Load movie
mdf_file='/Volumes/Analysis/movie-xml/RGB-8-1-0.48-11111.xml';
triggers=datarun.triggers;
frames=20000%floor(max(triggers)*120) % 10 minutes
[mov,height,width,duration,refresh] = get_movie(mdf_file, triggers,frames);
mov=(mov-0.5);

%% Load spikes
% Right now, uses all cells! TODO - careful!! ..

neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile('/Volumes/Analysis/2012-08-09-3/data002/data002.neurons');
neuronIDs = neuronFile.getIDList();
TTLtimes=double(neuronFile.getTTLTimes());

nCell=length(neuronIDs);
vision_spk_times=cell(nCell,1);

for iCell=1:nCell
vision_spk_times{iCell}= double(neuronFile.getSpikeTimes(neuronIDs(iCell)));
end

spkBinned = cell(nCell,1);
binFrames=cell(nCell,1);
FrameDivision=1;
binLen=(1/FrameDivision)*(1/120)*20000; % in Samples
noFramesUse=10000;
startFrame=100;

for iCell=1:nCell
    spkBinned{iCell}=zeros(noFramesUse*FrameDivision,1);    
    icnt=0;
    iCell
    for iFrame=startFrame:startFrame+noFramesUse
        for idivision=1:FrameDivision
            icnt=icnt+1;
            numTriggers=floor(iFrame/100)+1;
            sampleStart = TTLtimes(numTriggers) + (iFrame-(numTriggers-1)*100)*(20000/120)+(idivision-1)*(20000/(120*FrameDivision));
            sampleEnd = TTLtimes(numTriggers) + (iFrame-(numTriggers-1)*100)*(20000/120)+(idivision)*(20000/(120*FrameDivision));
            
            spkBinned{iCell}(icnt)=sum(double(vision_spk_times{iCell}>sampleStart) .* double(vision_spk_times{iCell}<=sampleEnd ));
            binFrames{iCell}(icnt)=iFrame;
        
        end
    end
end

%spkBinned
a=cell(length(matlab_cell_ids),1);
for icell=1:length(matlab_cell_ids)
a{icell}=spkBinned{matlab_cell_ids(icell)};
end
spkBinned=a;
%binFrames
binFrames=binFrames{1};

%% 


% Use Recursive Least Squares
% http://web.stanford.edu/class/ee263/lectures/06_ls-app.pdf - 6-21,
% Notation from there
delay=20; % Change this ?? according to movie length etc? 
mov_len=1000;
cellIds=matlab_cell_ids;
noCells=length(cellIds);

spk_coll=spkBinned;

filter_bank=cell(size(mov,1),size(mov,2),3);

for pixX=15%1:size(mov,1);
for pixY=12%1:size(mov,2);
for pixCol=1:1;
pixX
pixY
pixCol
mov_recons = squeeze(mov(pixX,pixY,pixCol,:));
q=zeros(delay*noCells+1,1);
P=zeros(delay*noCells+1,delay*noCells+1);
nSamples=mov_len*FrameDivision-delay;
A=zeros(nSamples,delay*noCells+1);
b=zeros(nSamples,1);
for iLen=1:nSamples%1800*120
    iFrame=binFrames(iLen)
    a=[1];% add a constant
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iLen:1:iLen+delay-1))];
    end
    a=full(a);
    y=mov_recons(iFrame);

    A(iLen,:)=a';
    b(iLen)=y;
    
    q=q+y*a;
    P=P+a*a';

end

 
filter=P\q;   
filter2=(A\(A'\(A'*b)));

filter_use=filter2;

% Make prediction
mov_pred=0*mov_recons;
recons_idx=mov_len*FrameDivision:size(mov,4)*FrameDivision-delay
for iLen=recons_idx %1800*120
    iFrame=binFrames(iLen)
    iFrame
    a=[1];
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iLen:1:iLen+delay-1))];
    end
    a=full(a);
    
   mov_pred(iLen)= filter_use'*a;
    
end

figure;
stairs(mov_recons(binFrames(recons_idx(1)+100:recons_idx(1)+200)),'b')
hold on
stairs(mov_pred(recons_idx(1)+100:recons_idx(1)+200),'r');
%hold on;
%stairs(double(mov_pred(recons_idx(1)+100:recons_idx(1)+200)>0)-0.5,'g');
xlim([1,100])
% correlation 

filter_bank{pixX,pixY,pixCol}.filter=filter_use;
filter_bank{pixX,pixY,pixCol}.correlation=corr(mov_pred(recons_idx),mov_recons(binFrames(recons_idx)));
corr(mov_pred(recons_idx),mov_recons(recons_idx))
corr(double(mov_pred(recons_idx)>0),mov_recons(recons_idx))
end
end
end



% save('/Volumes/Analysis/nishal/Reconstruct_2010-03-05-2_data018_1.mat','datarun');

% save('/Volumes/Analysis/nishal/Reconstruct_2010-03-05-2_data018_2.mat','filter_bank','mov','movie_spec','-v7.3');

