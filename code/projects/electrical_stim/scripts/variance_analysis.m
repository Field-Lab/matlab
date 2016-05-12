[rawData, amplitudes] = generateEiFromStimPattern('/Volumes/Analysis/2016-02-17-9/data003/', [451],'suppressPlots',true);

%firstArtifact = mean(rawData(:,:,:,1),1);

%subtractionMatrix = repmat(firstArtifact,[size(rawData,1) 1]);

%subtractionMatrix =repmat(subtractionMatrix,[1 1 1 size(rawData,4)]);

%modData = rawData - subtractionMatrix;
modData = rawData;

%var_modData = squeeze(var(modData,0,1));
var_modData = squeeze(var(modData));

%get amplitudes


 

%figure; 

for ii = [15]

   %hold on;  plot(squeeze(var_modData(15,:,ii)),'-','Color',C(ii,:),'LineWidth',2);
   figure; 
   %hold on;  plot(linspace(0,45,100),squeeze(var_modData(369,:,ii)),'LineWidth',2);
   plot(squeeze(var_modData(415,1:45,ii)),'LineWidth',2);

end
%{
for ii = [9 12 15]

   %hold on;  plot(squeeze(var_modData(15,:,ii)),'-','Color',C(ii,:),'LineWidth',2);
   figure; 
   %hold on;  plot(linspace(0,45,100),squeeze(var_modData(369,:,ii)),'LineWidth',2);
   hold on;  plot(squeeze(var_modData(369,1:45,ii)),'LineWidth',2);

end
%}

%{
%plot of maxiumum variance----------------------------------------
figure;
max_vals = squeeze(max(var_modData(336,1:45,:),[],2)); 
ampmat = [.27803 .30331 .35386 .37913 .42969 .45496 .50551 .55606 .60661 .68244 .73299 .80882 .88465 .98575 1.1043 1.2047 1.3051 1.4055 1.6063 1.7067 1.9075 2.1083 2.3091 2.5098 2.811 3.1122 3.4134 3.7146]; 
plot(ampmat, max_vals)
%}	

%summation analysis----------------------------------------
%{
figure
foo=squeeze(mean(modData(:, 369, :, 16),1))' 
plot(foo)
figure
foosum = 0;
for iii = 1:15; foosum = foosum + foo*normrnd(.75,.2); end
plot(foosum)
%}



%{
figure; 

for t = 1:size(var_modData,2)

     scatter(xc,yc,200,var_modData(:,t,28),'filled');colorbar; 

     caxis([0 500]); axis image; axis off; 

     title(num2str(t)); 

     pause(0.1);

end

mean_var = squeeze(mean(var_modData,2)); 

meanData = squeeze(mean(modData,1));

allAmps = squeeze(max(meanData,[],2) - min(meanData,[],2));

figure;

for a = 1:size(mean_var,2)

    

    subplot(2,1,1);

    cla;

    scatter(xc,yc,200,allAmps(:,a),'filled');colorbar;caxis([0 200]);

    hold on; scatter(xc([336]),yc([336]),200,'r','filled');

    axis image; axis off;

    subplot(2,1,2);

%    cla; 

%}
%{
    scatter(xc,yc,200,mean_var(:,a),'filled');colorbar;

    hold on; scatter(xc([336]),yc([336]),200,'r','filled');

    caxis([0 200]); axis image; axis off;

    title(num2str(a));

    pause(0.1);

end

elec = 322; 

elec = 15; 

elec = 314; 

elec = 327; 

figure; 

for ii = 23:38

    hold on; plot(squeeze(modData(:,elec,1:50,ii))','-','Color',C(ii,:));

    hold on; plot(squeeze(mean(modData(:,elec,1:50,ii)))','-','Color',0.5*C(ii,:),'LineWidth',2);

    title(num2str(ii)); 

    pause; 

end

figure; 

for ii = 1:38

   hold on;  plot(squeeze(var_modData(elec,:,ii)),'-','Color',C(ii,:),'LineWidth',2);

end

meanData = squeeze(mean(modData,1));

allAmps = squeeze(max(meanData,[],2) - min(meanData,[],2));
%}
