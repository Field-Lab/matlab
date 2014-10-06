function slope=przekladki_channels_scan(a,b);

dac=[-127:127];
slopea=ones(1,32)*nan;
slopeb=slopea;
for i=[1:32]
    s=a((i-1)*256+1:i*256)';
    [slope1,lin_err]=dac_lin5(dac,s);
    slopea(i)=slope1;    
    s=b((i-1)*256+1:i*256)';
    [slope2,lin_err]=dac_lin5(dac,s);    
    slopeb(i)=slope2;    
end

slope=slopeb./slopea;