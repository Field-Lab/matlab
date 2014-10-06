size(characteristics);

ch=zeros(n,length(frqs));
gen_ampl=250*0.001;
tlumik=800;

figure(1);
clf;
for i=1:64
    ch(:,:)=characteristics(i,:,:)/2/gen_ampl*tlumik; % 2 - bo odczyt roznicowy!
    subplot(8,8,i);
    %hold on;
    %for j=1:6
     semilogx(frqs,mean(ch));   
     grid on;
%    plot(frqs,ch(j,:));
%end
end    

figure(2);
clf;

for i=1:64
    ch(:,:)=characteristics(i,:,:)/2/gen_ampl*tlumik;
    %hold on;
    %for j=1:6
    semilogx(frqs,mean(ch),'k*-');   
    grid on;
    hold on;
%    plot(frqs,ch(j,:));
%end
end    

points=[1 4 6 8 10 13];  %numery czestotliwosci, dla ktorych robimy wykresy wzm. vs numer kanalu

for i=1:6
    s=zeros(64,6);
    s(:,:)=characteristics(:,:,points(i))/2/gen_ampl*tlumik;
    figure(3);
    subplot(3,2,i);
    plot(mean(s'),'k*-');
    %text(50,,'dghdfhg');
    grid on;
end