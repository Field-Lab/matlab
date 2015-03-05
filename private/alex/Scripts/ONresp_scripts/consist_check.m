function cons=consist_check(data)

tmp=triu(corr(data),1);
tmp=tmp(tmp~=0);
cons=nanmean(tmp);