% basis_maker
clear
load ps_bases.mat;
unit_basis = ps_bases{2}.basis;
bins = size(unit_basis,1);

%%
dense.index = [1:20];
dense.vector = zeros(bins,1);
dense.vector(1:2) = [1;1];
dense.rotation_factor = 1;
dense.rotation_0 = 1; 

sparse.index = [1:19];
sparse.vector = zeros(bins,1);
sparse.vector(1:10) = ones(10,1);
sparse.rotation_factor = 5;
sparse.rotation_0      = 20;

new_basis = [];

for i_portion = 1:2
    if i_portion == 1, pars = dense; end
    if i_portion == 2, pars = sparse; end
    
    new_part = [];
    for i_ind = pars.index
        new_vec = circshift(pars.vector, ...
            pars.rotation_0 + (i_ind-1)*pars.rotation_factor);
        new_part = [new_part new_vec];
    end
    
    new_basis = [new_basis new_part];
end

imagesc(new_basis);

ps_bases{3}.basis = new_basis;
ps_bases{3}.note  = 'steps_overlap_40';

save('ps_bases.mat', 'ps_bases')
%%

dense.index = [1:20];
dense.vector = zeros(bins,1);
dense.vector(1) = 1;
dense.rotation_factor = 1;
dense.rotation_0 = 1; 

sparse.index = [1:45];
sparse.vector = zeros(bins,1);
sparse.vector(1:10) = ones(10,1);
sparse.rotation_factor = 2;
sparse.rotation_0      = 20;

new_basis = [];

for i_portion = 1:2
    if i_portion == 1, pars = dense; end
    if i_portion == 2, pars = sparse; end
    
    new_part = [];
    for i_ind = pars.index
        new_vec = circshift(pars.vector, ...
            pars.rotation_0 + (i_ind-1)*pars.rotation_factor);
        new_part = [new_part new_vec];
    end
    
    new_basis = [new_basis new_part];
end

imagesc(new_basis);

ps_bases{4}.basis = new_basis;
ps_bases{4}.note  = 'smooth_steps_overlap_65';

save('ps_bases.mat', 'ps_bases')

%%
dense.index = [1:20];
dense.vector = zeros(bins,1);
dense.vector(1) = 1;
dense.rotation_factor = 1;
dense.rotation_0 = 1; 

sparse.index = [1:18];
sparse.vector = zeros(bins,1);
sparse.vector(1:5)  = (linspace(0,.8,5))'; 
sparse.vector(6:15) = ones(10,1);
sparse.vector(16:20) = (linspace(.8,0,5))';
sparse.rotation_factor = 5;
sparse.rotation_0      = 15;

new_basis = [];

for i_portion = 1:2
    if i_portion == 1, pars = dense; end
    if i_portion == 2, pars = sparse; end
    
    new_part = [];
    for i_ind = pars.index
        new_vec = circshift(pars.vector, ...
            pars.rotation_0 + (i_ind-1)*pars.rotation_factor);
        new_part = [new_part new_vec];
    end
    
    new_basis = [new_basis new_part];
end

imagesc(new_basis);

ps_bases{5}.basis = new_basis;
ps_bases{5}.note  = 'rampedsteps_overlap_39';

save('ps_bases.mat', 'ps_bases')

