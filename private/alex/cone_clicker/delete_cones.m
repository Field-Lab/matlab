function delete_cones

global cone_list h i

myCone=round(ginput(1));
a=triu(squareform(pdist([cone_list; myCone])));
a(a==0)=max(a(:));
[r,c]=find(a==min(a(:)));
myCone=min(r,c);

if (ishandle(h(myCone)))
    delete(h(myCone));
    h(myCone)=[];
    i=i-1;
    cone_list(myCone,:)=[];
end
 




