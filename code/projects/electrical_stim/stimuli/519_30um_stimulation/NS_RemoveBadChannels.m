function GoodChannels=NS_RemoveBadChannels(Channels,BadChannels);

GoodChannels=[];
for i=1:length(Channels)
    active=1;
    for j=1:length(BadChannels)
        if Channels(i)==BadChannels(j)
            active=0;
        end
    end
    if active==1
        GoodChannels=[GoodChannels Channels(i)];
    end
end