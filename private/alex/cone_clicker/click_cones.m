function click_cones

global cone_list h myPicture
global i

if isempty(i);i=1;end

hold on
xsize=get(gca,'Xlim');
ysize=get(gca,'Ylim');
flag=1;
while flag
    myCone=round(ginput(1));
    if myCone(1)>=xsize(1) && myCone(1)<=xsize(2) && myCone(2)>=ysize(1) && myCone(2)<=ysize(2)
        %optimize location
        mySub=myPicture(myCone(2)-5:myCone(2)+5, myCone(1)-5:myCone(1)+5);
        [a,b]=find(mySub==max(mySub(:)),1);
        myCone(2)=myCone(2)+a-6;
        myCone(1)=myCone(1)+b-6;
        cone_list=[cone_list; myCone];
        h(i)=plot(myCone(1),myCone(2),'rx');
        i=i+1;
    else
        flag=0;
    end
end

