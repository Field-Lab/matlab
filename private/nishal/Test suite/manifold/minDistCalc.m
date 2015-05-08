function E_d_row = minDistCalc(E,Edge01,icity)     

 
    idx=1:size(E,1);
E_d = E(icity,:);
E_d(E_d==0 & Edge01(icity,:) == logical(0))=Inf;

un_visited =idx(E_d==Inf);
un_visited = un_visited(un_visited~=icity);

while(length(un_visited)~=0)
icnt=0;
unvisited_d=[];
for tovisit=un_visited
    
    icnt=icnt+1;
    neighbors = idx(Edge01(:,tovisit)~=0);
    
    d_min=Inf;
   
    for ineighbor =neighbors
       
        d=E_d(ineighbor)+ E(ineighbor,tovisit);
        if(d<d_min)
        d_min=d;
        min_neigh=ineighbor;
        end
    end
    
    unvisited_d(icnt) = d_min;

end
    length(unvisited_d);
[v,i]=min(unvisited_d);
add_city = un_visited(i);
E_d(add_city)=v;
un_visited =idx(E_d==Inf);
un_visited = un_visited(un_visited~=icity);
length(un_visited);

E_d_row=E_d;
end