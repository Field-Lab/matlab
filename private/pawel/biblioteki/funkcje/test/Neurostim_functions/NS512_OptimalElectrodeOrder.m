function OrderFor512=NS512_OptimalElectrodeOrder();

%We picture the array as having 16 rows and 32 columns (we ignore that it
%is hexagonal and not rectangular). The first index is row and the second
%is column.

OrderFor16(1:16,1)=[1 1 3 3 2 2 4 4 1 1 3 3 2 2 4 4]; %4 x 4 array
OrderFor16(1:16,2)=[1 3 1 3 1 3 1 3 2 4 2 4 2 4 2 4];
OrderFor256=zeros(256,2);
%Now order for 256 electrodes (16 x 16 array)
for electrode=1:16 %number of 16-electrode cluster...
    for cluster=1:16 %number of electrode within the culster...
        index=(electrode-1)*16+cluster;
        OrderFor256(index,1)=(OrderFor16(cluster,1)-1)*4+OrderFor16(electrode,1);
        OrderFor256(index,2)=(OrderFor16(cluster,2)-1)*4+OrderFor16(electrode,2);
    end
end

OrderFor512=zeros(512,2);
OrderFor512(1:2:511,1)=OrderFor256(:,1);
OrderFor512(2:2:512,1)=OrderFor256(:,1)+16;
OrderFor512(1:2:511,2)=OrderFor256(:,2);
OrderFor512(2:2:512,2)=OrderFor256(:,2);