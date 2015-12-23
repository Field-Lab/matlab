function [v,c]=ei_axon_speed(datarun,cell_identifier,varargin)
%finding centroid in time
%greschner

% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('plot', 0);% 
    p.addParamValue('threshold', 2);% 
    p.addParamValue('cutoff', .3);%cutoff thr for values % smaller than max 
    p.addParamValue('nr_average', 1);%
    p.addParamValue('min_nr_electrodes', 4);%
    p.addParamValue('axon_electrodes', []);%
    p.addParamValue('fit_method', 'single');%'single' fit or 'multi'ple fits or 'accelerate'ing 
    
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;



if params.plot
    subplot(1,2,1)
end

if isempty(params.axon_electrodes)
    res=ei_get_axon(datarun,cell_identifier,'plot',params.plot,'threshold',params.threshold);
else
    res=params.axon_electrodes;
end


v=NaN; c=NaN;
if length(res)>=params.min_nr_electrodes    
    index=get_cell_indices(datarun,cell_identifier);

    ei=datarun.ei.eis{index}(res,:);
    po=datarun.ei.position(res,:);

    tei=zeros(size(ei));
    for i=1:size(ei,1)
        [t1,t2]=sort(ei(i,:));
        %tei(i,t2(end-params.nr_average+1:end))=t1(end-params.nr_average+1:end);
        tei(i,t2(1:params.nr_average))=-t1(1:params.nr_average);
    end

    tei=tei/max(abs(tei(:)));
    tei(tei<params.cutoff)=0;

    c=zeros(2,size(ei,2));
    for i=1:size(ei,2)
        ei_frame=tei(:,i);
        if any(ei_frame)
            c(1,i)=sum(ei_frame.*po(:,1))/sum(ei_frame);
            c(2,i)=sum(ei_frame.*po(:,2))/sum(ei_frame);
        else
            c(1,i)=NaN;
            c(2,i)=NaN;
        end
    end

    t=find(~isnan(c(1,:)));
    tp=t-t(1);
    tc=c(:,t);
    
    if params.plot
        hold on
        plot(c(1,:),c(2,:),'.-')
        hold off
    end

    
    
    %test move forward
    for i=2:length(tp)
        theta(i)=cart2pol(tc(1,i)-tc(1,i-1),tc(2,i)-tc(2,i-1));
    end
    tt=[];
    for i=3:length(tp)
        if abs(theta(i)-theta(i-1))>3/4*pi
            disp(sprintf('centroid does move back - delete last %d frames',length(tp)-i+1));
            tt=[tt i];
        end
    end
    if ~isempty(tt)
        tp=tp(1:tt(1)-1);
        tc=tc(:,1:tt(1)-1);
    end

    if params.plot
        hold on
        plot(tc(1,:),tc(2,:),'g.-')
        hold off
    end
 
    tc=tc/1000; %transform mm
    tp=tp/20; %transform ms
    

    

    if isequal(params.fit_method, 'single')
        if length(tp)>=params.min_nr_electrodes
            %m1=regress(di',tp');
            %dtp=tp(2:end)-tp(1:end-1);
            %m3=robust_mean(ddi./dtp);
            
            for i=1:length(tc)-1
                di(i)=sqrt((tc(1,i+1)-tc(1,i))^2+(tc(2,i+1)-tc(2,i))^2);
            end 
            di=[0 cumsum(di)];
            
            v=robustfit(tp,di,'bisquare',4.685,'off');

            if params.plot
                subplot(1,2,2)
                plot(tp,di,'.');
                hold on
                plot(tp,tp*v,'-r');
                hold off

                textloc(sprintf('%.3f',v),'NorthWest');
            end
        end
    end


    if isequal(params.fit_method, 'multi')

        if length(tp)>=params.min_nr_electrodes
            vv=[];
            for i=1:length(tc)-params.min_nr_electrodes
                r=zeros(params.min_nr_electrodes,2); 
                for ii=1:params.min_nr_electrodes
                    r(ii,1)=sqrt((tc(1,ii+i-1)-tc(1,i))^2+(tc(2,ii+i-1)-tc(2,i))^2);
                    r(ii,2)=tp(ii+i-1)-tp(i);
                end               
                m=robustfit(r(:,2),r(:,1),'bisquare',4.685,'off');
                vv=[vv m];

                if params.plot
                    col=colormap('lines');
                    subplot(1,2,2)
                    plot(r(:,2),r(:,1),'.','color',col(i,:));
                    hold on
                    plot(r(:,2),r(:,2)*m,'-','color',col(i,:));
                end
            end 

            v=mean(vv);

            if params.plot
                textloc(sprintf('%.3f',v),'NorthWest');
                hold off
            end
        end
    end

    
    if isequal(params.fit_method, 'accelerate')

        if length(tp)>=params.min_nr_electrodes
            for i=1:length(tc)-1
                r(i,1)=sqrt((tc(1,i+1)-tc(1,i))^2+(tc(2,i+1)-tc(2,i))^2);
                r(i,2)=tp(i+1)-tp(i);
                r(i,3)=sqrt((tc(1,1)-tc(1,i))^2+(tc(2,1)-tc(2,i))^2);
                
                vv(i)=r(i,1)/r(i,2);
            end 

            v=mean(vv);
     
            if params.plot
                subplot(1,2,2)
                plot(r(:,3),vv,'.-');
                hold on
                textloc(sprintf('%.2f',mean(vv)),'NorthWest');
                hold off
            end

        end
    end
    
    
    
end