%X=NS512_OptimalElectrodeOrder();
%X=NS512_OptimalElectrodeSequence();
%X=NS512_RowAndColumnForElectrode(500,PatternsForSet)

Disconnected=[1 130 259 260 389 390 519];
GoodChannels=NS_RemoveBadChannels(Electrodes,Disconnected);

ElectrodeOrder=zeros(1,512);
ElectrodeOrder(1:128)=GoodChannels(1:4:512);
ElectrodeOrder(129:256)=GoodChannels(2:4:512);
ElectrodeOrder(257:384)=GoodChannels(3:4:512);
ElectrodeOrder(385:512)=GoodChannels(4:4:512);

X0=electrode_positions(519);
X=X0(ElectrodeOrder,:);

a=X(:,1);
b=X(:,2);

figure(3);
axis([0 33 0 17]);
clf;
for i=1:519
    plot(a(1:i,1),b(1:i,1),'bd');
    axis([-400 400 -400 400]);
    hold on;
    pause(0.1);
    %axis([0 33 0 17]);
end