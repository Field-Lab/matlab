



nthres=5;
flags=zeros(nNeurons,1);

while(prod(flags)==0)
    flags=zeros(nNeurons,1);
    conddelet=zeros(E,J); 
    clear caux
   for n=1:nNeurons
        caux=[];
       for j=1:J
       if(nansum(spikes{n}(j,:)==I(j))
           caux=[caux j];


if(~isempty(caux))
    cdel(ne)=caux(1);
    suma=nansum(nansum(spikes{ne}(1:caux(1)-1,:)'));
    
    if(suma<nthres)
conddelet(ne,cdel(ne))=-1;
   
else
    flags(ne)=1;
end
else
    flags(ne)=1;
   end
   end
   