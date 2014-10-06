function y=plot_deconv(animal,type,ktore_spiki,system,filtr_odwr,dl,figura);
%Rysuje kilka spikow - na jednym wykresie: z "war data" i zdekonwoluowanych, dla okreslonego zwierzaka i typu spikow.
%animal: 1 - malpa, 2 - guinepig
%type: 1 - cell body, 2 - axonal, ew. kiedys dalsze
%ktore spiki - numery spikow - nie przegiac!
%system: 1 - linux, 2 - windows
%filtr_odwr - y2 zwracana przez inveyefilter_hayes
%dl: parametr podawany takze funkcji inveyefilter_hayes (aby skorygowac przesuniecie)
%figura - na ktorym wykresie rysowac.

if animal==1
	[p,filenames,spikes,detect_param]=spiki_malpa(type,ktore_spiki,system);
else
	[p,filenames,spikes,detect_param]=spiki_gunpg(type,ktore_spiki,system);
end
	
%dl=200;
%[y1,y2]=inveyefilter_hayes(dl,58,40,3300,20000);

margins=[15 50];
t=[0:margins(1)+margins(2)]/20;

ls=length(ktore_spiki);
sqls=sqrt(ls);
if floor(sqls)*floor(sqls)==ls
	sb1=floor(sqls)
	sb2=floor(sqls);
else
	if floor(sqls)*ceil(sqls)>=ls
		sb1=floor(sqls)
		sb2=ceil(sqls);
	else
		sb1=ceil(sqls);
		sb2=ceil(sqls);
	end
end
	

figure(figura);
for j=1:length(filenames)
    if find(spikes(:,1)==j)'
        s=importdata(filenames{j})';
        s=s-mean(s);
    
        sy=conv(s,filtr_odwr);
        sy=sy(dl+1:length(s)+dl);
        sy=sy-mean(sy);

        positions=find(spikes(:,1)==j)' %numery dobrych spikow w tym pliku    
        %coordinates(positions)=y(1,spikes(positions,2)); 

        for i=1:length(positions)
	    subplot(sb1,sb2,positions(i));
	    strt=spikes(positions(i),2)-margins(1)
	    stp=spikes(positions(i),2)+margins(2)
	    plot(t,s(strt:stp),'b-',t,sy(strt:stp),'r-');
	    grid on
        end
    end    
end


