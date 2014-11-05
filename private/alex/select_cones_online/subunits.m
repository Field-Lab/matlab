function subunits(datarun,path2load, date)


global cellType myCells cones stim

datarun = load_neurons(datarun);

datarun.names.rrs_movie_path=[path2load,'.movie'];

% leave only selected cells

for i=1:length(myCells)
    myIndices(i)=find(datarun.cell_ids==myCells(i));
end
datarun.cell_ids=datarun.cell_ids(myIndices);
datarun.cell_types{cellType}.cell_ids=myCells;
datarun.channels=datarun.channels(myIndices);
datarun.spikes=datarun.spikes(myIndices);

% faking cones
allcones=cell2mat(cones');
datarun.cones.centers=allcones;
datarun.cones.types=repmat('U',size(allcones));
datarun.cones.weights=zeros(length(allcones),length(myCells));
cnt=1;
for i=1:length(myCells)
    for j=1:size(cones{i},1)
       datarun.cones.weights(cnt,i)=1;
       cnt=cnt+1;
    end
end 
 
Wc = sparse(prod(datarun.stimulus.field_height*datarun.stimulus.field_width)*3,length(allcones)); % the hardcoded 3 makes this an image


kern_size=size(stim.cone_template,1)-1;
tmp_coord=stim.coord-kern_size/2;


for i=1:length(allcones);
    
    bw_kern=zeros(datarun.stimulus.field_height,datarun.stimulus.field_width);
    bw_kern(tmp_coord(i,2):tmp_coord(i,2)+kern_size,...
        tmp_coord(i,1):tmp_coord(i,1)+kern_size)=stim.cone_template(:,:,stim.templateNumber(i));
    % reshape to a single column
    bw_kern = reshape(bw_kern,[],1);
    
    % get rgb values
    rgb = [1 1 1];
    % multiply out
    rgb_kern = bw_kern * rgb;
    % reshape
    rgb_kern = reshape(rgb_kern,[],1);
    % normalize sum
    rgb_kern = rgb_kern/sum(rgb_kern);
    % put into W
    Wc(:,i) = (rgb_kern);
end

% load('/Users/alexth/Desktop/myWc.mat');
% Wc=myWc;
% datarun.cones.centers=coords;

tmp=regexp(datarun.names.rrs_prefix,'/');

save([datarun.names.rrs_prefix(1:tmp(end)),'Wc.mat'],'Wc');
save([datarun.names.rrs_prefix(1:tmp(end)),'cone_info.mat'],'cones','myCells', 'stim');

datarun = conepreprocess_wnm(datarun, 'conepath',datarun.names.rrs_prefix(1:tmp(end)));

conepreprocess_save(datarun, 'cell_types',datarun.cell_ids,'cone_data_ind', 'data', 'date', date);



% %%
% stixel=datarun.stimulus.stixel_height;
% kern=imresize(cone_templ,1/stixel);
% kern_size=floor(size(kern)/2);
% tmp_coord=round(stim.coord/stixel);
% 
% for i=1:length(allcones);
%     
%     bw_kern=zeros(datarun.stimulus.field_height,datarun.stimulus.field_width);
%     bw_kern(tmp_coord(i,2)-kern_size(2):tmp_coord(i,2)+kern_size(2),...
%         tmp_coord(i,1)-kern_size(1):tmp_coord(i,1)+kern_size(1))=kern;
%     % reshape to a single column
%     bw_kern = reshape(bw_kern,[],1);
%     
%     % get rgb values
%     rgb = [1 1 1];
%     % multiply out
%     rgb_kern = bw_kern * rgb;
%     % reshape
%     rgb_kern = reshape(rgb_kern,[],1);
%     % normalize sum
%     rgb_kern = rgb_kern/sum(rgb_kern);
%     % put into W
%     Wc(:,i) = (rgb_kern);
% end
% % figure
% % imagesc(reshape(full(Wc(1:90000,i)),300,300))


