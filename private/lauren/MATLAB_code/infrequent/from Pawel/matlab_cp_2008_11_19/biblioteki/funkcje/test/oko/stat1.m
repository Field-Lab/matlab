function [srednia,odchylenie]=channel_stat(dane);
r=size(dane);

for i=1:r(1)
   srednia(i)=mean(dane(i,:));
   odchylenie=std(dane(i,:));
end
