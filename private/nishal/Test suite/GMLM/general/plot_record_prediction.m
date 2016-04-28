function h3=plot_record_prediction(spkCondColl,spkCondCollGLM,varargin)

p = inputParser;
% specify list of optional parameters
PS_str.plot=0;
p.addParamValue('PSTH_struct', PS_str);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% assign variables from parsing.
PSTH_str = params.PSTH_struct;




col='rkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkrkr';
nmov = length(spkCondColl);

dt=1/120;

xx=cell(nmov,1);
yy=cell(nmov,1);
for imov=1:nmov

%plotraster(x,fittedGLM,'labels',true,'raster_length',24)
[xx{imov},yy{imov}]=plotSpikeRaster(spkCondCollGLM{imov}>0,'PlotType','vertline');
xx{imov}=xx{imov}*dt;

end


scrsz = get(groot,'ScreenSize');
h3=[];
% h3= figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);

jump=0;
icnt=0;
for imov=nmov:-1:1
    
    % Simulated
   icnt=icnt+1;
    plot(xx{imov},yy{imov}+jump,col(icnt));%,'Marker','o','MarkerSize',1);
    hold on;
    if(PSTH_str.plot~=0)
        PSTH=PSTH_str.PSTH_pred_log{imov};PSTH = PSTH-mean(PSTH(:));PSTH =( PSTH/max(abs(PSTH(:)))*max(yy{imov}(:))/2)+max(yy{imov}(:))/2;
        plot([0:dt:(length(PSTH)-1)*dt],jump+PSTH,col(icnt+1),'LineWidth',2);
    end
    jump=jump+30.5;%max(yy{imov});
    

    % Recorded 
    icnt=icnt+1;
    plot(spkCondColl(imov).xPoints/20000,spkCondColl(imov).yPoints+jump,col(icnt));%,'Marker','o','MarkerSize',1);
    hold on
    if(PSTH_str.plot~=0)
        PSTH=PSTH_str.PSTH_rec_log{imov};PSTH = PSTH-mean(PSTH(:));PSTH =( PSTH/max(abs(PSTH(:)))*max(spkCondColl(imov).yPoints(:))/2)+max(spkCondColl(imov).yPoints(:))/2;
        plot([0:dt:(length(PSTH)-1)*dt],jump+PSTH,col(icnt+1),'LineWidth',2);
    end
    jump = jump + max(spkCondColl(imov).yPoints);



end
ylim([0,jump]);
xlim([0,max(xx{imov})]);
set(gca,'YTick',[]);
end