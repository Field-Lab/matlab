function e=getNeighbors(e1,deg)


for k=1:deg
    aux=e1;
    for j=1:length(e1)
        aux=union(getCluster512(e1(j)),aux);
    end
    e1=aux;
end
e=e1;