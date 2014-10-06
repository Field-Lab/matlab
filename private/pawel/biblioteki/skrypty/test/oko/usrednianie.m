start=1;
dlugosc=10000000;
%r=readconv('Data001conv',206,65,[8 9 10 18 19],[(start) (start+dlugosc-1)]);

kroki=500;
okno=10000;

r=zeros(5,kroki);

for i=1:kroki
   c=readconv('Data000conv',206,65,[8 9 10 18 19],[(i*okno+1) ((i+1)*okno)]);
   d=sum(c')';
   r(1:5,i)=d;
end


figure(1);
subplot(3,3,2);
plot(r(1,:));
%axis([0 dlugosc -400 400]);
subplot(3,3,4);
plot(r(2,:));
%axis([0 dlugosc -400 400]);
subplot(3,3,5);
plot(r(3,:));
%axis([0 dlugosc -400 400]);
subplot(3,3,6);
plot(r(4,:));
%axis([0 dlugosc -400 400]);
subplot(3,3,8);
plot(r(5,:));
%axis([0 dlugosc -400 400]);