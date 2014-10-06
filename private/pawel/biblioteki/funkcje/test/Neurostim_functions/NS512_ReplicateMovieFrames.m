function OutputMovie=NS512_ReplicateMovieFrames(Movie,N);

l=length(Movie);

for i=1:l
    index=(i-1)*N;
    OutputMovie(index+1:index+N)=Movie(i);
end