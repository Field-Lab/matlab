function y=readxls(filename,start,length);

y=textread(filename,'',length,'headerlines',start);