function file_time
%displays timestamp of calling functions (rehash problem by calls from terminal) 
%greschner

t1=dbstack;
t2=which(t1(end).file);
t3=dir(t2);
t4=24*60*60*(now-datenum(t3.date));

display(sprintf('\n%s was saved %s - %.0f sec ago\n',t1(end).file,t3.date,t4));