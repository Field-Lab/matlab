function frame=get_frame(n, seed, dim)

defaultStream = RandStream.getGlobalStream;
defaultStream.State = seed.State;

frame=zeros(numel(n),dim);

if numel(n)==1
    for i=1:n
        frame=randi([0 1],1,dim);
    end
else
    for i=1:n(1)-1
        randi([0 1],1,dim);
    end
    j=1;
    for i=n(1):n(2)
        frame(j,:)=randi([0 1],1,dim);
        j=j+1;
    end
end
end