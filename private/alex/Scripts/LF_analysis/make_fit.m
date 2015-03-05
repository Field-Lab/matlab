function make_fit(currentLFs, start_trial,trialsPerScreen,fitParams,currentUnit)

global fitResults

end_trial=min(start_trial+trialsPerScreen-1,size(currentLFs,2));

resFit=zeros(3,end_trial-start_trial+1);

cnt=1;

for trial=start_trial:end_trial
    if fitParams.SpikeCount(trial)>5
        borders=fitParams.leftBorder(trial):fitParams.rightBorder(trial);
        
        currentFit=currentLFs(borders,trial);        
        currentFit=currentFit*fitParams.currentPolarity;
        
        fit_res=fit(borders',currentFit,'gauss1',...
            'Lower',fitParams.limits(1,1:3),'Upper',fitParams.limits(2,1:3),...
            'Startpoint',[fitParams.amp(trial), fitParams.ind(trial), 30]);
        resFit(1,cnt)=fitParams.currentPolarity*fit_res.a1;
        resFit(2,cnt)=fit_res.b1;
        resFit(3,cnt)=fit_res.c1;
    end
    cnt=cnt+1;
end

fitResults(:,start_trial:end_trial,currentUnit)=resFit;
makePlots(currentLFs,trialsPerScreen,start_trial,resFit,0);

display('DONE')

