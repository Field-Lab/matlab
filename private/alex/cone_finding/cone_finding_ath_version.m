
% datarunA=load_sta(datarunA,'load_sta','all');

cell_types={1,2,3,4};
cell_ids=datarunA.cell_ids(get_cell_indices(datarunA, cell_types));



datarunA.names.nickname = '';
datarunA.piece.rig = 'A';
datarunA.piece.optical_path_direction = 'below';
datarunA.piece.display = 'crt1';
extra_dirname_info = 'RGB-1-8';

datarunA = set_polarities(datarunA);


% BW or RGB stimulus?
independent = strcmpi(datarunA.stimulus.independent, 't');
field_width = datarunA.stimulus.field_width;
field_height = datarunA.stimulus.field_height;

% robust_std_method is 1 to match old implementation.  Set to 3,5 for some speedup.
disp('Loading sta summaries and static NLs...')
datarunA = get_sta_summaries(datarunA, cell_types, ...
    'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct( ...
        'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
        'thresh',5,'robust_std_method',1));
    
% calculate static nonlinearities
datarunA = load_java_movie(datarunA, '/Volumes/Analysis/deprecated/movie-xml2/RGB-1-8-0.48-11111-320x320.xml');
start_time=0;
datarunA = get_snls(datarunA, datarunA.cell_ids(get_cell_indices(datarunA, cell_types)),'frames',-2:0,'start_time',start_time,'stimuli',10000,'new',true);

fit_b=zeros(1,length(cell_ids));
cell_constants=fit_b;
n_spikes=fit_b;
cell_cr=fit_b;
for i=1:length(cell_ids)
    fit_b(i) = datarunA.stas.snls{datarunA.cell_ids==cell_ids(i)}.fit_params.b;
    cell_constants(i)=exp(fit_b(i));    
    n_spikes(i) = length(datarunA.spikes{datarunA.cell_ids==cell_ids(i)});
    cell_cr(i)=(n_spikes(i) / cell_constants(i));
end



resize_factor=3;
pad=10*resize_factor;

sta_all=zeros(320*resize_factor,320*resize_factor,length(cell_ids));
kk=zeros(pad*2+1,pad*2+1,length(cell_ids));
cnt=1;


for i=1:length(cell_ids)
    tmp=sum(datarunA.stas.stas{datarunA.cell_ids==cell_ids(i)}(:,:,:,6),3);
    tmp=imresize(tmp,resize_factor);
    [a,b]=max(abs(tmp(:)));
    
    sta_all(:,:,i)=tmp*sign(tmp(b))*cell_cr(i);
    
    tmp=tmp/a;
    [x,y]=ind2sub(size(tmp),b);
    if y>pad && x>pad && x<(320*resize_factor-pad) && y<(320*resize_factor-pad)
        kk(:,:,cnt)=tmp(x-pad:x+pad,y-pad:y+pad)*sign(tmp(b));
        cnt=cnt+1;
    end

end

kk=kk(:,:,1:cnt-1);

figure
subplot(2,1,1)
tmp=mean(kk,3);
imagesc(tmp/max(tmp(:)))
subplot(2,1,2)
tmp=std(kk,0,3);
imagesc(tmp)

tmp=mean(kk,3);
tmp=tmp-robust_mean(tmp(:));
normCone=tmp/max(tmp(:));
figure
imagesc(normCone)


myPrior=[4.0,4.25]; % L-M cones prior

myCone=normCone(24:38,24:38);

figure;imagesc(myCone)
myCones=repmat(myCone,1,1,size(sta_all,3));
p=zeros(size(sta_all));
tic
for i=(size(myCone,1)+1):size(sta_all,1)-(size(myCone,1)+1)
    for j=(size(myCone,1)+1):size(sta_all,2)-(size(myCone,1)+1)
        tmp=sta_all((1:size(myCone,1))+i,(1:size(myCone,2))+j,:).*myCones;
        p(i+ceil(size(myCone,1)/2),j+ceil(size(myCone,1)/2),:)=sum(reshape(tmp,numel(myCone),153));
    end
end
toc


figure
colormap gray
subplot(2,1,1)
tmp=reshape(p,960*960,153);
tmp=max(tmp');
tmp=reshape(tmp,960,960);
imagesc(tmp/max(tmp(:)));
% axis([440 540 100 200 ])
robust_mean(tmp(:)/max(tmp(:)))
subplot(2,1,2)
colormap gray
tmp=reshape(sta_all,960*960,153);
tmp=max(tmp');
tmp=reshape(tmp,960,960);
imagesc(tmp/max(tmp(:)));
% axis([440 540 100 200 ])
robust_mean(tmp(:)/max(tmp(:)))

figure
colormap gray
tmp=reshape(p,960*960,153);
tmp=sort(tmp','descend');
a=tmp>0.55;
tmp=tmp(1,:)+tmp(2,:).*a(2,:);
tmp=reshape(tmp,960,960);
imagesc(tmp/max(tmp(:)));


figure
colormap gray
imagesc(fliplr(dllA))

figure
imagesc(p(:,:,1))




figure
imagesc(dllA)

myCone=normCone(24:38,24:38);
myCone=imresize(myCone,3);
figure
imagesc(myCone)

newCones=FastPeakFind(dllA,0.1,myCone,1,1);
figure
colormap gray
imagesc(dllA); hold on
plot(newCones(1:2:end),newCones(2:2:end),'r+')
hold on
plot(size(dllA,1)-conesA(:,1),conesA(:,2),'yx')

plot(newCones(1:2:end)+1510+9,newCones(2:2:end)+435,'yx')
plot(size(dllA,1)-conesA(:,1)+1510+9,conesA(:,2)+435,'+','color',[0.2 0.7 0.5])


m=pdist2([size(dllA,1)-conesA(:,1),conesA(:,2)],[newCones(1:2:end),newCones(2:2:end)]);

nnd=min(m');
figure
hist(nnd, 0:0.25:35)


tmp=reshape(p,960*960,153);
tmp=max(tmp');
tmp=reshape(tmp,960,960);
tmp=tmp/max(tmp(:));
tmp=fliplr(tmp);
tmp=imresize(tmp,3);
figure
imagesc(tmp);


myCone=imresize(myCone,3);
figure
imagesc(myCone)


newCones=FastPeakFind(tmp,0.1,myCone);
figure
colormap gray
imagesc(tmp); hold on
plot(newCones(1:2:end),newCones(2:2:end),'r+')
hold on
plot(size(dllA,1)-conesA(:,1),conesA(:,2),'yx')
m=pdist2([size(dllA,1)-conesA(:,1),conesA(:,2)],[newCones(1:2:end),newCones(2:2:end)]);

nnd=min(m');
figure
hist(nnd, 0:0.25:35)
