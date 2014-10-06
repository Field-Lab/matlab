function [ ChannelsPlot_mod ] = getAdjacent_modified( ChannelsPlot )
% do_eliminacji - Macierz zawieraj�ca te elektrody, kt�re eliminujemy z analizy(elektrody uszkodzone) w celu poprawienia wynik�w 
do_eliminacji=[9,25,57,4,31];
new=double(ChannelsPlot);
for i=1:max(size(do_eliminacji))
if find(new==do_eliminacji(i))
    pos=find(new==do_eliminacji(i));
    new=cat(2,new(1:pos-1),new(pos+1:end));
end
ChannelsPlot_mod=new;
end

