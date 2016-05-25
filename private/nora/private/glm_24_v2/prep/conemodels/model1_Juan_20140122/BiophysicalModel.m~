%% Created by Juan Angueyra (Jan_2014)
path('/netapp/snle/lab/temp_matlabcodeAH/fromJuan_v1_20140122/BiophysicalModel',path) 
%% Fit of biophysical model to dim flash response
load(sprintf('/Users/akheitman/Dropbox/ConeModelsAndGLM/BiophysicalModel/%s.mat','DimFlash_092413Fc12'))
dimflashfit=-riekemodelfx2(DimFlash.coeffs,DimFlash.timeaxis);

f1=getfigH(1);
set(get(f1,'XLabel'),'String','Time (s)')
set(get(f1,'YLabel'),'String','i (pA)')

lH=line(DimFlash.timeaxis,DimFlash.dimflashfit,'Parent',f1);
set(lH,'Color','r','LineStyle','-','Marker','.','DisplayName','fit');
lH=line(DimFlash.timeaxis,DimFlash.dimflash,'Parent',f1);
set(lH,'Color','k','LineStyle','-','Marker','.','DisplayName','dimflash');
lH=line(DimFlash.timeaxis,dimflashfit,'Parent',f1);
set(lH,'Color','b','LineStyle','none','Marker','o','DisplayName','fit');

%% Fit to extra-simplified "EyeMovement" stimulus adding extra calcium feedback
load(sprintf('/Users/akheitman/Dropbox/ConeModelsAndGLM/BiophysicalModel/%s.mat','EyeMovements_092413Fc12'))

fit=-riekemodelstim_slowadapt([DimFlash.coeffs(5) EyeMovements.slowcacoeff],EyeMovements.timeaxis,EyeMovements.stim,DimFlash.coeffs(1:7));
fit_nocaslow=-riekemodelstim_slowadapt([DimFlash.coeffs(5) DimFlash.coeffs(8)],EyeMovements.timeaxis,EyeMovements.stim,DimFlash.coeffs(1:7));


f2=getfigH(2);
set(get(f2,'XLabel'),'String','Time (s)')
set(get(f2,'YLabel'),'String','i (pA)')

lH=line(EyeMovements.timeaxis,EyeMovements.data,'Parent',f2);
set(lH,'Color','k','DisplayName','response');
lH=line(EyeMovements.timeaxis,EyeMovements.fit_nocaslow,'Parent',f2);
set(lH,'Color','b','LineStyle','-','DisplayName','No_CaSlow');
lH=line(EyeMovements.timeaxis,EyeMovements.fit,'Parent',f2);
set(lH,'Color','r','LineStyle','-','DisplayName','CaSlow');
lH=line(EyeMovements.timeaxis,fit_nocaslow,'Parent',f2);
set(lH,'Color','c','Marker','none','LineStyle',':','DisplayName','No_CaSlow');
lH=line(EyeMovements.timeaxis,fit,'Parent',f2);
set(lH,'Color','m','Marker','none','LineStyle',':','DisplayName','No_CaSlow');




stimnormed = (EyeMovements.stim - min(EyeMovements.stim)) / max(EyeMovements.stim - min(EyeMovements.stim));
fitnormed  =  (fit - min(fit)) / max(fit - min(fit));

t = EyeMovements.timeaxis;
figure;
subplot(3,1,1); plot( t ,  EyeMovements.stim, 'b');
subplot(3,1,2); plot( t , fit, 'r');
subplot(3,1,3); hold on; plot(t ,stimnormed,'b'); plot(t,fitnormed, 'r')


%%