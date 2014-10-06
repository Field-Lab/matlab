function fill_up

global myCells cones daratun

all_cones=cell2mat(cones');

figure

plot(all_cones(:,1),all_cones(:,2),'*r')
axis ij
hold on


a=convhull(all_cones(:,1),all_cones(:,2));


IN = inpolygon(datarun.cones.centers(:,1),datarun.cones.centers(:,2),...
    datarun.cones.centers(all_COI(a),1),datarun.cones.centers(all_COI(a),2));

%plot results
plot(all_cones(a,1),all_cones(a,2),'r-')

,...
    datarun.cones.centers(IN>0,1),datarun.cones.centers(IN>0,2),'g*',...
    datarun.cones.centers(all_COI,1),datarun.cones.centers(all_COI,2),'m*')
%     plot(datarun.cones.centers(IN>0,1),datarun.cones.centers(IN>0,2),'g+',...
