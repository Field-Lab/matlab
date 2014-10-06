function Subsubplot(WholeAreaCoordinates,RelativeSubplotCoordinates);
NewLeft=WholeAreaCoordinates(1)+RelativeSubplotCoordinates(1)*WholeAreaCoordinates(3);
NewBottom=WholeAreaCoordinates(2)+RelativeSubplotCoordinates(2)*WholeAreaCoordinates(4);
NewWidth=WholeAreaCoordinates(3)*RelativeSubplotCoordinates(3);
NewHeight=WholeAreaCoordinates(4)*RelativeSubplotCoordinates(4);
subplot('position',[NewLeft NewBottom NewWidth NewHeight]);