function TimingsForSpikesMinimumsCorrected=NS512_CorrectSpikesTimings(TimingsForSpikesMinimums,SpikesWithinTimeBrackets);

TimingsForSpikesMinimumsCorrected=TimingsForSpikesMinimums;
g2=min(TimingsForSpikesMinimums);
if g2>4
    TimingsForSpikesMinimumsCorrected=TimingsForSpikesMinimums-4;
else
    TimingsForSpikesMinimumsCorrected=TimingsForSpikesMinimums-g2+1;
end

%{
g1=find(SpikesWithinTimeBrackets==1)     
TimingsForSpikesMinimumsCorrected=TimingsForSpikesMinimums
if length(g1)>0 % jesli sa jakies spiki z niewielkim jitterem
    g2=min(TimingsForSpikesMinimums(g1))
    if g2>4
        TimingsForSpikesMinimumsCorrected(g1)=TimingsForSpikesMinimums(g1)-4;
    else
        TimingsForSpikesMinimumsCorrected(g1)=TimingsForSpikesMinimums(g1)-g2+1;
    end
end    
%}