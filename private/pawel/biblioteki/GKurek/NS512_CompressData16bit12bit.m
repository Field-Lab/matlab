function out=NS512_CompressData16bit12bit(OutputData);

s1 = uint16(OutputData(2:2:end,:)+2048);
s2 = uint16(OutputData(3:2:end,:)+2048);
    
out=zeros(770,size(OutputData,2),'uint16');

b1 = bitshift(s1,-4);       
b2 = bitshift(bitand(s1,15),4) + bitshift(s2,-8);       
b3=bitand(s2,255);

out(3:3:end,:) = b1;
out(4:3:end,:) = b2;
out(5:3:end,:) = b3;   