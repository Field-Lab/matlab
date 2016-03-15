function y= lastConnectedComponent(x)
%Gonzalo Mena, 3/2016
cont=1;

aux=x(1);
for i=2:length(x)
    if(x(i)==(aux(end)+1))
        aux(end)=x(i);
    else
        cont=cont+1;
        aux(cont)=x(i);
    end
end
     y=aux;   