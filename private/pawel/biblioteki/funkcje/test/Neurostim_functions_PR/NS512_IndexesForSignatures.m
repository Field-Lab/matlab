function Indexes=NS512_IndexesForSignatures(SpikesTimings,TimeCorelSpikes);

g1=find(TimeCorelSpikes==1);         
if ~isempty(length(g1)) % jesli sa jakies spiki z niewielkim jitterem
    g2=min(SpikesTimings(g1));
    if g2>4
        SpikesTimings(g1)=SpikesTimings(g1)-4;
        else
            SpikesTimings(g1)=SpikesTimings(g1)-g2+1;
        end
end   
Indexes=g1;