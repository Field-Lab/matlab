function ElectrodeOrder=NS512_519OptimalOrder;
ec=electrode_positions(519);

x=round(ec(:,1));
y=round(ec(:,2));

Electrodes=[];

for xc=[min(x):30:max(x)]
    xc
    el0=find(x==xc);
    yc=y(el0)
    [af,i]=sort(yc);
    Electrodes=[Electrodes' el0(i)']';
end

Disconnected=[1 130 259 260 389 390 519];
GoodChannels=NS_RemoveBadChannels(Electrodes,Disconnected);

ElectrodeOrder=zeros(1,512);
ElectrodeOrder(1:128)=GoodChannels(1:4:512);
ElectrodeOrder(129:256)=GoodChannels(2:4:512);
ElectrodeOrder(257:384)=GoodChannels(3:4:512);
ElectrodeOrder(385:512)=GoodChannels(4:4:512);