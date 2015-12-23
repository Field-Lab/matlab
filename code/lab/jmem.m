function jmem

max_mem = java.lang.Runtime.getRuntime.maxMemory / 2^20;
tot_mem = java.lang.Runtime.getRuntime.totalMemory / 2^20;
free_mem = java.lang.Runtime.getRuntime.freeMemory / 2^20;

fprintf('%0.1f of %0.1f MB used (%d%%), %0.1f MB max\n',...
    tot_mem-free_mem,tot_mem,round(100*(tot_mem-free_mem)/tot_mem),max_mem)
