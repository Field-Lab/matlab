function ConvertedData=ChannelDataToColors2(ChannelData,X,Y);
%ChannelData - one value per each 
%X,Y - coordinates for each electrode

ConvertedData=zeros(64,16);

Indexes=NS512_RowAndColumnForElectrode(500,[1:512]);

for i=1:512
    Xind=round([X(i)/30-0.5 X(i)/30+0.5])+33;
    Yind=max(Indexes(:,2))-Indexes(i,2)+1;
    ConvertedData(Xind,Yind)=ChannelData(i);
end    