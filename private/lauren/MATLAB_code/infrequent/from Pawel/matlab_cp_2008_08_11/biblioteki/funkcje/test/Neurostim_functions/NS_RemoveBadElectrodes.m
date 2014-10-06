function GoodChannels=NS_RemoveBadChannels(Channels,BadChannels);

GoodElectrodes=[];
for i=1:length(Electrodes)
    active=1;
    for j=1:length(BadElectrodes)
        if Channels(i)==BadChannels(j)
            active=0;
        end
    end
    if active==1
        GoodChannels=[GoodChannels Channels(i)];
    end
end