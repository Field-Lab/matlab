function [data_len,op1,op2] = collect_ends(data)
data_len=[];op1=[];op2=[];
for itime = 1:length(data)
    data_len(itime) = data(itime).useMin;
    fval_test = data(itime).fval_test; fval_test=gather(fval_test);
    d_true = (data(itime).d_true); d_true = gather(d_true);
    op1(itime) =fval_test(end);
    op2(itime) =d_true(end);
end