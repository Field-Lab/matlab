function update_cone_info

global cones myCells mean_nnd

all_cones=[];
for i=1:length(cones)
    all_cones=[all_cones; cones{i}];
end
all_cones=squareform(pdist([all_cones(:,1),all_cones(:,2)]));
all_cones(all_cones==0)=10000;
all_cones=min(all_cones);
mean_nnd=mean(all_cones);

hInfo=uicontrol('style','text', 'Units', 'Normalized','position',[0.32 0.05 0.16 0.08],...
    'string',{[int2str(length(cell2mat(cones'))),' cones in ', int2str(length(myCells)), ' cells'],'',...
    ['Mean nnd: ', num2str(mean_nnd,2),' +/- ', num2str(std(all_cones),2)]},...
    'fontsize',16, 'fontweight', 'bold');
