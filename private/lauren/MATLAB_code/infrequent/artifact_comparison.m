close all

vectorLength = 161;

figure
subplot(6,1,1)
title('Artifact Comparison for Stimulaton with Electrode 14, 2.2 µA')
hold on
for j = 1:7
	plot(1 + (j-1)*vectorLength : j*vectorLength, data009failureAvg(j,:),'r-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data029failureAvg(j,:),'m-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data036failureAvg(j,:),'b-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data043failureAvg(j,:),'k-')
%      for i = 1:size(data009failuresToSave, 1)
%          plot(1 + (j-1)*vectorLength : j*vectorLength, squeeze(data009failuresToSave(i,j,:))-3000,'r-')
%      end
end
hold off



subplot(6,1,2)
title('Artifact Comparison for Stimulaton with Electrode 14, 2.6 µA')

hold on
for j = 1:7
	plot(1 + (j-1)*vectorLength : j*vectorLength, data010failureAvg(j,:),'r-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data030failureAvg(j,:),'m-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data037failureAvg(j,:),'b-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data044failureAvg(j,:),'k-')
end
hold off


subplot(6,1,3)
title('Artifact Comparison for Stimulaton with Electrode 19, 1.2 µA')

hold on
for j = 1:4
	plot(1 + (j-1)*vectorLength : j*vectorLength, data014failureAvg(j,:),'r-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data031failureAvg(j,:),'m-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data038failureAvg(j,:),'b-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data045failureAvg(j,:),'k-')
end
hold off


subplot(6,1,4)
title('Artifact Comparison for Stimulaton with Electrode 19, 1.5 µA')

hold on
for j = 1:4
	plot(1 + (j-1)*vectorLength : j*vectorLength, data015failureAvg(j,:),'r-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data032failureAvg(j,:),'m-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data039failureAvg(j,:),'b-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data046failureAvg(j,:),'k-')
end
hold off



subplot(6,1,5)
title('Artifact Comparison for Stimulaton with Electrode 54, 2.6 µA')

hold on
for j = 1:7
	plot(1 + (j-1)*vectorLength : j*vectorLength, data025failureAvg(j,:),'r-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data033failureAvg(j,:),'m-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data040failureAvg(j,:),'b-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data047failureAvg(j,:),'k-')
end
hold off



subplot(6,1,6)
title('Artifact Comparison for Stimulaton with Electrode 54, 3.0 µA')

hold on
for j = 1:7
	plot(1 + (j-1)*vectorLength : j*vectorLength, data026failureAvg(j,:),'r-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data034failureAvg(j,:),'m-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data041failureAvg(j,:),'b-')
	plot(1 + (j-1)*vectorLength : j*vectorLength, data048failureAvg(j,:),'k-')
end
hold off

legend('0 nM TTX', '50 nM TTX', '100 nM TTX', '1 µM TTX')