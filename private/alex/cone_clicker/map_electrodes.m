function map_electrodes

persistent  hh
if ~isempty(hh) && ishandle(hh(1))
    for i=1:length(hh)
        if ishandle(hh(i))
            delete(hh(i))
        end
    end    
end 
hh=[];
    
myEl=round(ginput(2));
if myEl(1,2)>myEl(2,2)
    myEl=[myEl(2,:); myEl(1,:)];
end

maxDist=pdist(myEl);
singleDist=maxDist/26;
borderDistance=0.46*maxDist;


beginsX=(myEl(1,1)-borderDistance):borderDistance/12:(myEl(1,1)+borderDistance);
beginsY=[(myEl(1,2)+0.5*borderDistance):(-0.5*borderDistance)/12:myEl(1,2) (myEl(1,2)+0.5*borderDistance/12):(0.5*borderDistance)/12:(myEl(1,2)+0.5*borderDistance) ];
endsY=[(myEl(2,2)-0.5*borderDistance):0.5*borderDistance/12:myEl(2,2) (myEl(2,2)-0.5*borderDistance/12):-0.5*borderDistance/12:(myEl(2,2)-0.5*borderDistance) ];
hold on

for i=1:25
    myCol=beginsY(i):singleDist:(endsY(i)+singleDist/2);
    hh(i)=plot(repmat(beginsX(i),1,length(myCol)),myCol,'xy');
end
 

 