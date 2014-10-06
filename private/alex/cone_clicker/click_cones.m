function click_cones(sbh,hclick)

global cone_list myPicture

hold on
xsize=get(gca,'Xlim');
ysize=get(gca,'Ylim');
flag=1;
while flag
    myCone=round(ginput(1));
    if myCone(1)>=xsize(1) && myCone(1)<=xsize(2) && myCone(2)>=ysize(1) && myCone(2)<=ysize(2) 
        cone_list=[cone_list; myCone];
        plot(myCone(1),myCone(2),'rx')
    else
        flag=0;
    end
end

