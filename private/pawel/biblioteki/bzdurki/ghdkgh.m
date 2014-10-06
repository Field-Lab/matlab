figure(10)
clf
for n=47:512
%n=256;
n
clf
subplot(1,2,1);
h=plot(RawData(n,:));
set(h,'Color','b');
hold on
h=plot(ArtifactData(n,:));
set(h,'Color','r');
grid on
axis([0 10000 -400 -100])

subplot(1,2,2);
plot(RawData(n,:)-ArtifactData(n,:));
axis([0 10000 -100 100])
grid on
pause(1)
end