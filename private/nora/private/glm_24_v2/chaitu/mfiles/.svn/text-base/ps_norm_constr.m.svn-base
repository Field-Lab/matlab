function [c ceq GC GCeq] = ps_norm_constr(x,basepars,stimpars,trainpars)

if (~isfield(basepars,'psnormconstr_p'))
    basepars.psnormconstr_p = 2; %p  in Lp norm constraint
end
p = basepars.psnormconstr_p;

ceq = 0;
GCeq = zeros(size(x));

c = zeros(basepars.Nneurons,1);
GC = zeros(size(x,1),basepars.Nneurons);

for j=1:basepars.Nneurons
    ps_idx = get_pars_idx(basepars,j,size(trainpars.D,2),'ps');
    c(j) = 1/p * (sum(abs(x(ps_idx)).^p) - basepars.psnormconstr^p);
    GC(ps_idx,j) = abs(x(ps_idx)).^(p-1).*sign(x(ps_idx));
end
