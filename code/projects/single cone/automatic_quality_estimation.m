myData=[];
myRuns={}; 
j=1;k=0;
for i=1:length(WNruns)
    a=WNruns{i};
    if length(a)>3 && strcmp(a(1:4),'data')
        myData=[myData i];
        if ~isempty(dataset1{i})   
            j=1;
            k=k+1;
            myRuns{k,1}=dataset1{i};
            myRuns{k,2}{j}=a;            
        else
            j=j+1;
            myRuns{k,2}{j}=a;           
        end
    end
end


for myDataruns=7:length(myRuns)
    piece=myRuns{myDataruns,1};
    for myData00=1:length(myRuns{myDataruns,2})
        run = myRuns{myDataruns,2}{myData00};
        
        % define data path
        myPath=['/Volumes/Analysis/' piece '/' run ];
        if ~exist(myPath,'dir')
            myPath=['/Volumes/Analysis/' piece '/streamed/' run ];
        elseif exist([myPath,'/',run],'dir')
            myPath=[myPath,'/',run];
        end
        fid=fopen('/Users/alexth/Desktop/dataset_quality.txt','a');
        fprintf(fid,'\n%s\t%s\t',piece,run);
        
        if exist([myPath,'/',run,'.sta'],'file')
            datarunA = load_data([myPath,'/',run]);
            opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',false,'load_time_courses',true);
            datarunA=load_data(datarunA,opt);
            tic
            datarunA = load_sta(datarunA,'load_sta','all');
            toc
            datarunA=get_sta_summaries(datarunA,'all');
            
            
            fram=datarunA.stas.depth-1;
            figure
            cnt=1;
            crap=length(datarunA.cell_types);

            for j=[1:5 crap]
                myCells=datarunA.cell_types{1,j}.cell_ids;
                myCone=zeros(41);
                for i=1:length(myCells)
                    tmp=find(datarunA.cell_ids==myCells(i));
                    if datarunA.stimulus.independent=='t'
                        if j~=5
                            a=sum(datarunA.stas.stas{tmp}(:,:,:,fram),3)/3;
                        else
                            a=squeeze(datarunA.stas.stas{tmp}(:,:,3,fram));
                        end
                    else
                        a=squeeze(datarunA.stas.stas{tmp}(:,:,1,fram));
                    end
                    minn=min(a(:));
                    maxx=max(a(:));
                    if abs(minn)>maxx %OFF cell
                        a=-a;
                    end
                    
                    [val,loc]=max(a(:));
                    [x,y]=find(a==val,'1');
                    if x>20 & x<datarunA.stimulus.field_height-20 & y>20 & y<datarunA.stimulus.field_width-20
                        myCone=myCone+a(x-20:x+20,y-20:y+20);
                    end
                end
                subplot(2,3,cnt)
                imagesc(myCone)
                tmp=myCone(21,21)/(sum((sum(abs(myCone(20:22,20:22)))))-myCone(21,21));
                tmp2=((sum((sum(abs(myCone(20:22,20:22)))))-myCone(21,21))/8)/(sum(sum(abs(myCone(1:10,1:10))))/100);
                title({[datarunA.cell_types{1,j}.name, '  n=',int2str(length(myCells)),],[num2str(tmp),'  snr=',num2str(tmp2)]})
                fprintf(fid,'%d\t%.2f\t%.2f\t',length(myCells), tmp, tmp2);
                cnt=cnt+1;
            end
            saveas(gcf,['/Users/alexth/Desktop/single_cone_datasets_quality/',piece,'_',run,'.png'])
            close gcf
            clear datarunA myPath 
        else
            fprintf(fid,'%s','data or sta not found');
        end
        fclose(fid);
    end
end    


myVars=who;

for j=1:length(myVars)
    for i=1:106
        if isempty(eval([myVars{j},'{i}']))||eval([myVars{j},'{i}(1)'])=='d'||eval([myVars{j},'{i}(1)'])=='u'||eval([myVars{j},'{i}(1)'])=='p'||eval([myVars{j},'{i}(1)'])=='t'
            eval([myVars{j},'{i}=NaN']);
        elseif ~isnan(eval([myVars{j},'{i}']))
            eval([myVars{j},'{i}=str2num(eval([myVars{j},''{i}'']))'])
        end
    end
end
for j=1:length(myVars)
    eval([myVars{j},'=cell2mat(',myVars{j},');']);
end

for i=2:length(tempecc)
    if isnan(tempecc(i))
        tempecc(i)=tempecc(i-1);
    end
end


plot(tempecc,SBCnumber,'*')
hold on
plot(tempecc,ONmidgetnumber,'*r')
plot(tempecc,OFFmidgetnumber,'*g')
plot(tempecc,ONparasolnumber,'*m')
plot(tempecc,OFFparasolnumber,'*c')

figure
plot(tempecc,SBCSNRbiggerbetter,'*')
hold on
plot(tempecc,ONmidgetSNRbiggerbetter,'*r')
plot(tempecc,OFFmidgetSNRbiggerbetter,'*g')
plot(tempecc,ONparasolSNRbiggerbetter,'*m')
plot(tempecc,OFFparasolSNRbiggerbetter,'*c')


