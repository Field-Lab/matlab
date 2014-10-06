s3=s2+25;
znak=sign(s3);

s=size(s3);
znaki=zeros(s);
for i=1:s(1)
    z=znak(i,:);
    znaki(i,:)=z;
    figure(11)
    subplot(5,4,i)
    plot(z,'bd-')
    figure(12);
    subplot(5,4,i);
    plot(diff(z),'bd-');
    figure(14)
    k=find(diff(z)<0);
    %subplot(5,4,i);
    h=plot(k,[i],'bd');
    set(h,'MarkerSize',2);
    hold on;
end
    

d=diff(znak);
find(d==-1);