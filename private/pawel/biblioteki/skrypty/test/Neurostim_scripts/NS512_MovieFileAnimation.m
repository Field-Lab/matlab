Chunks=MovieChunksFileCombined;

NumberOfMovies=Chunks(1);

figure(101);
clf;
index=2;
for i=1:NumberOfMovies
    length=Chunks(index)
    %if (find(Movies==i) | Movies==0)
        data=Chunks(index+7:index+length);
        sd=size(data);
        patterns=data(2:3:sd(2))
        h=NS512_electrode_map_animation2(patterns);
    %end
    index=index+length+1;
end