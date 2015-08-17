  
for nSU=1:7
load(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data008_gmlm_fig/quad_gmlm_3736_su_%d_data.mat',nSU));
    close all
    
h=figure;
            for icone = 1:nc
                for jcone=1:nc
                    iconeidx =cones(icone); jconeidx=cones(jcone); pair = [iconeidx;jconeidx];
                    if(su_log(icone,jcone)>0 & icone~=jcone)
                        plot(datarun.cones.centers(pair,1),-datarun.cones.centers(pair,2),'LineWidth',su_log(icone,jcone));
                        hold on;
                    end
                end
            end
            
            for icone = 1:nc
                for jcone=1:nc
                    iconeidx =cones(icone); jconeidx=cones(jcone); pair = [iconeidx;jconeidx];
                    if(su_log(icone,jcone)>0 & icone~=jcone)
                        x1 = mean(datarun.cones.centers(pair,1));
                        y1 = mean(-datarun.cones.centers(pair,2));
                        text(x1,y1,sprintf('%d',su_log(icone,jcone)));
                        hold on;
                    end
                end
            end
            
            for icone=1  :nc
                scatter(datarun.cones.centers(cones(icone),1) ,-datarun.cones.centers(cones(icone),2),abs(sta(icone,end-1))*2000,'filled','r');
                hold on;
            end
            xlim([min(datarun.cones.centers(cones,1))-2,max(datarun.cones.centers(cones,1))+2]);
            
            ylim([min(-datarun.cones.centers(cones,2))-2,max(-datarun.cones.centers(cones,2))+2])
            title(sprintf('no. of SU: %d',nSU));
            
            set(gca,'xTick',[]);
            set(gca,'yTick',[]);
           hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data008_gmlm_fig/quad_gmlm_%d_su_%d.eps',cellID,nSU));
end

 %%
  x=repmat(1:7,[7,1]);
  y=repmat([1:7]',[1,7]);
  
 mask = logical(x>y)
for nSU=1:7
load(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data008_gmlm_fig/quad_gmlm_3736_su_%d_data.mat',nSU));
    close all
ss_l = su_log(mask);
h=figure;
hist(ss_l,25);
xlim([0,55]) 
set(gca,'yTick',[]);
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data008_gmlm_fig/quad_gmlm_hist_3736_su_%d.eps',nSU));
pause
end