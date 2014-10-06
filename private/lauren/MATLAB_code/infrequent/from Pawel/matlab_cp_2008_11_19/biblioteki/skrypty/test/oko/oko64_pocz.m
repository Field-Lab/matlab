cd /home/pawel/pliki/oko/oko_64/pomiary_pocz/

list=['tek00000.dat'; 'tek00001.dat'; 'tek00006.dat'; 'tek00007.dat'];
list=['tek00010.dat'; 'tek00011.dat'; 'tek00012.dat'; 'tek00013.dat'; 'tek00014.dat'; 'tek00015.dat'; 'tek00016.dat'; 'tek00017.dat'; 'tek00018.dat';];
colors=['b'; 'r'; 'g'; 'k'];
a=size(list);

figure(6);
figure(7);
clf(6);
clf(7);
for i=1:4
    a0=importdata(list(i,:));
    m=mean(a0((length(a0)-1000):length(a0),1));
    s=a0(5:length(a0),1)';
    b(i,:)=s;
    %figure(i+1);
    %subplot(2,2,i);
    kolor=i-floor((i-1)/length(colors))*length(colors);
    plot(s,colors(kolor));
    hold on;
end
grid on;