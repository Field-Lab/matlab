
%% structure of program
% load datarun and sta (all classified cells? or specified cell type?)
% explore possibility of automatic classification at this point.
% plot rf fits for specified cell type
% click on cells (within the rf fit) to show sta (zoomed in) in a subplot
% (possibly several at once? several suplots if space allows)
% for each sta, look for cones automatically (use max and subtract some
% area around). If clicked again in same cell, remove it

% Parameters fields (editable in GUI): area around the cone
% (in stixels), max distance from the nearest neighbor (in stixels, should
% be recalculated at each iteration: if a new cone is suddenly much farther
% away from the rest, stop; allow for manual adjustment of the parameter.
% Perform checking for max allowed distance starting from first 2 cones -
% if it's not SBC, it should not be too big, otherwise crap).

% printed info: N of cones, mean distance, show average cone template?

% button 'adjust positions': click it and get best cone position according
% to gaussian template (or average data template if desired). Then, draw a
% region around each cone's RF to show which area will be stimulated.

% button 'verify other cell types': look for cells with RF fits in the same
% area and check if their cone positions correspond. Adjust if desired?

% button 'add/correct manually': if clicked, add cone manually; if cone
% already exists close by, correct position; if cone already exists in
% exact location, delete it

% Show combined plot of cone locations from a group of cells in a separate
% suplot/main window where rf fits were shown?

% button 'save results': saves maps for stimulation.



%%


 datarun = load_data(fullfile(server_path(), '2013-08-19-2/streamed/data001/data001'));
 
 
 
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);


datarun = load_sta(datarun,'load_sta','all');
field_width = datarun.stimulus.field_width;
field_height = datarun.stimulus.field_height;
bord=10;
cnt=0;
coneTempl=0;
for i=1:length(datarun.cell_ids)
    tmp=double(squeeze(datarun.stas.stas{i}));
    if abs(min(tmp(:)))>max(tmp(:))  % OFF cell, invert polarity
        tmp=-tmp;
    end
    myMax=max(tmp(:));
    myCoord=find(tmp==myMax,1);
    [x y z]=ind2sub(size(tmp),myCoord);
    if  x<(bord+1)||x>(field_width-bord)||y<(bord+1)||y>(field_height-bord)
        cnt=cnt+1;
    else
        coneTempl=coneTempl+tmp(x-bord:x+bord,y-bord:y+bord,z);
    end
end
    
figure
imagesc(coneTempl)

myCells=get_cell_indices(datarun, {4});
for i=1:length(myCells)
    tmp=double(squeeze(datarun.stas.stas{myCells(i)}));
    if abs(min(tmp(:)))>max(tmp(:))  % OFF cell, invert polarity
        tmp=-tmp;
    end
    
    figure
    colormap gray
    imagesc(tmp(:,:,4))
    
    myX=[];myY=[];
    for j=1:10
        myMax=max(tmp(:));
        myCoord=find(tmp==myMax,1);
        [x y z]=ind2sub(size(tmp),myCoord);
        myX=[myX,x];
        myY=[myY,y];
        tmp(x-2:x+2,y-2:y+2,:)=0;
    end
    hold on
    plot(myY,myX,'xr','markersize',20)
    
end
 
 
 %%
 movie_spec='/Volumes/Analysis/deprecated/movie-xml2/BW-3-3-0.48-11111-320x320.xml';

datarun.names.nickname = '';
datarun.piece.rig = 'A';
datarun.piece.optical_path_direction = 'below';
datarun.piece.display = 'crt1';
extra_dirname_info = 'RGB-1-6';


datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);


% BW or RGB stimulus?
independent = strcmpi(datarun.stimulus.independent, 't');
field_width = datarun.stimulus.field_width;
field_height = datarun.stimulus.field_height;

info(datarun);

cell_types = {1,2,3,4,5};


myCells=[47 486 556 601 617 878 887 1503 1847 1970 2132 3646 3904 5387 7656];
cellN = get_cell_indices(datarun, myCells);
datarun = load_sta(datarun,'load_sta',myCells);

% cheap autocorrelation
coneTempl=0;
bord=10;
for i=1:length(myCells)
    tmp=datarun.stas.stas{cellN(i)};
    a=tmp(:,:,2,5);
    [x,y]=find(a==max(a(:)));
    if x<(bord+1)||x>(320-bord)||y<(bord+1)||y>(320-bord)
        disp('skip')
        i
    else
        
        a=tmp(x-bord:x+bord,y-bord:y+bord,:,:);
        coneTempl=coneTempl+a;
    end
end

coneTempl=coneTempl/(length(myCells)-2);

a=min(coneTempl(:))*2;
coneTempl=coneTempl/a+0.5;
figure
for i=1:6
    subplot(2,3,i)
    imagesc(coneTempl(:,:,2,i))
end

figure
colormap gray
imagesc(coneTempl(:,:,2,5))


tmp=datarun.stas.stas{cellN(1)};

a=min(tmp(:))*2;
tmp=tmp/a+0.5;

figure
for i=1:6
    subplot(2,3,i)
    imagesc(tmp(:,:,:,i))
end
k=tmp(:,:,2,5);
[b,ic]=sort(k(:));

figure
hist(k(:),100)

robust_std(k(:))*3+robust_mean(k(:))

figure
plot(a(ic(end-200:end),:,2)')

figure
a=reshape(tmp,320*320,3,6);
a=permute(a,[1,3,2]);
plot(a(:,:,1)','r')
hold on
plot(a(:,:,2),'g')
plot(a(:,:,3),'b')

t=[];
for i=cellN
    tmp=datarun.stas.stas{i};
    a=min(tmp(:));
    tmp=tmp/a;
    k=tmp(:,:,:,5);
    k=sum(k,3);
    k=k(:);
    t=[t k];
end                                                             
figure
hist(t(:),100)

border=20;

figure
colormap gray
cnt=1;
for i=cellN
    tmp=datarun.stas.stas{i};
    k=-tmp(:,:,:,5);%+tmp(:,:,:,3)+tmp(:,:,:,2);
    k=k(:,:,2);
%     k=sum(k,3);
    k=k/max(k(:));
    min(k(:));
    t=3*robust_std(k(:));
    k(k<t)=0;
    [a,b]=find(k==1);
    subplot(4,4,cnt)
    imagesc(k)
    title(int2str(myCells(cnt)))
    axis([b-border b+border a-border a+border])
    cnt=cnt+1;
end




