%%%%% 2012-09-27-3

head_cellID = [6858 544  5192 31 1];
Neighbors_cellID = struct( 'head_ID',num2cell(head_cellID),'OFF_Par_neighbor_id',[],'ON_Par_neighbor_id',[]);
%%% BOTH NEIGHORS NOW INCLUDE ONE ODDBALL CELL THAT'S HEALTHY BUT SHOULDN'T
%%% BE THER CUZ NOT CLOSE
Neighbors_cellID(1).OFF_Par_neighbor_id = [6841 6856 6991 7006 7096]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(1).ON_Par_neighbor_id  = [3 35 139 6754 6995 7217] ; 
Neighbors_cellID(1).control_id    = []; %[ 2027 3736];

Neighbors_cellID(2).OFF_Par_neighbor_id = [196 301 511 616 781 797]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(2).ON_Par_neighbor_id  = [246 273 437 680 784 845 981] ; 
Neighbors_cellID(2).control_id    = []; %[ 2027 3736];

Neighbors_cellID(3).OFF_Par_neighbor_id = [4877 4952 5117 5266 5417 5566 5611 5671 ]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(3).ON_Par_neighbor_id  = [5087 5179 5419 5614 5627] ; 
Neighbors_cellID(3).control_id    = []; %[ 2027 3736];

Neighbors_cellID(4).OFF_Par_neighbor_id = [167 7082 7096 7351 7441 7591]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(4).ON_Par_neighbor_id  = [35 91 273 7217 7456 7562] ; 
Neighbors_cellID(4).control_id    = []; %[ 2027 3736];

Neighbors_cellID(5).OFF_Par_neighbor_id = [301 5566 5671 6377 6616 6841]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(5).ON_Par_neighbor_id  =  [139 246 5627 6243 6618 6858]; 
Neighbors_cellID(5).control_id    = []; %[ 2027 3736];





%%%%% 2012-08-09-3

head_cellID = [1426 3226]
Neighbors_cellID(1).OFF_Par_neighbor_id = [497 1292 ]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(1).ON_Par_neighbor_id  = [841 1352 1502 1772] ; 
Neighbors_cellID(1).control_id    = []; %[ 2027 3736];

Neighbors_cellID(2).OFF_Par_neighbor_id = [497 1292 ]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(2).ON_Par_neighbor_id  = [1352 1426 121 650] ; 
Neighbors_cellID(2).control_id    = []; %[ 2027 3736];

Neighbors_cellID(3).OFF_Par_neighbor_id = [2431 3361 3676]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(3).ON_Par_neighbor_id  = [2313 2926 3152 3647 3799] ; 
Neighbors_cellID(3).control_id    = []; %[ 2027 3736];

Neighbors_cellID(4).OFF_Par_neighbor_id = [49 167 556]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(4).ON_Par_neighbor_id  = [302 631 7681] ; 
Neighbors_cellID(4).control_id    = []; %[ 2027 3736];


%%%%;]
%2012-08-09-3
head_cellID = [7066 1786 1276 1471];
%Neighbors_cellID = struct( 'head_ID',num2cell(head_cellID),'OFF_Par_neighbor_id',[],'ON_Par_neighbor_id',[]);
%%% BOTH NEIGHORS NOW INCLUDE ONE ODDBALL CELL THAT'S HEALTHY BUT SHOULDN'T
%%% BE THER CUZ NOT CLOSE
Neighbors_cellID(1).OFF_Par_neighbor_id = [467 1471 5086 6256 7066 7292]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(1).ON_Par_neighbor_id  = [586 6257 6721 7068] ; 
Neighbors_cellID(1).control_id    = []; %[ 2027 3736];
 
Neighbors_cellID(2).OFF_Par_neighbor_id = [1292 1321 1471 1771 2506]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(2).ON_Par_neighbor_id  = [1276 1352 1772 2101] ; 
Neighbors_cellID(2).control_id    = []; %[ 2027 3736];
 
Neighbors_cellID(3).OFF_Par_neighbor_id = [1321 1471 1786 2506 1292 1771]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(3).ON_Par_neighbor_id  = [586 650 1205 1352 2101 2641] ; 
Neighbors_cellID(3).control_id    = []; %[ 2027 3736];
 
 
 
Neighbors_cellID(4).OFF_Par_neighbor_id = [467 1321 1786 2506 5086 7066]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(4).ON_Par_neighbor_id  = [586 1205 1276 4216 6721] ; 
Neighbors_cellID(4).control_id    = []; %[ 2027 3736];


%%%%% 2012-08-21-1

head_cellID = [ 3707 6783 6107 47];
Neighbors_cellID(1).OFF_Par_neighbor_id = [3153 3377 3994 4371 4593 4596]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(1).ON_Par_neighbor_id  = [106 3196 3931] ; 
Neighbors_cellID(1).control_id    = []; %[ 2027 3736];

Neighbors_cellID(2).OFF_Par_neighbor_id = [5716 5899 6637 7370]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(2).ON_Par_neighbor_id  = [6107 6470 6813] ; 
Neighbors_cellID(2).control_id    = []; %[ 2027 3736];

Neighbors_cellID(3).OFF_Par_neighbor_id = [5716 5899 6637 7370]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(3).ON_Par_neighbor_id  = [5311 5868 6061 6470 6813] ; 
Neighbors_cellID(3).control_id    = []; %[ 2027 3736];

Neighbors_cellID(4).OFF_Par_neighbor_id = [49 167 556]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(4).ON_Par_neighbor_id  = [302 631 7681] ; 
Neighbors_cellID(4).control_id    = []; %[ 2027 3736];





head_cellID = [1471];
Neighbors_cellID(1).OFF_Par_neighbor_id = [467 1321 1786 2506 5086 7066]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(1).ON_Par_neighbor_id  = [586 1205 1276 4216 6721] ; 
Neighbors_cellID(1).control_id    = []; %[ 2027 3736];

%{
Neighbors_cellID(2).OFF_Par_neighbor_id = [5716 5899 6637 7370]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(2).ON_Par_neighbor_id  = [6107 6470 6813] ; 
Neighbors_cellID(2).control_id    = []; %[ 2027 3736];

Neighbors_cellID(3).OFF_Par_neighbor_id = [5716 5899 6637 7370]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(3).ON_Par_neighbor_id  = [5311 5868 6061 6470 6813] ; 
Neighbors_cellID(3).control_id    = []; %[ 2027 3736];

Neighbors_cellID(4).OFF_Par_neighbor_id = [49 167 556]; %  ON [3 139 246 5627 6243 6618 6858]  %% DON'T INCLUDE THE NEIGHBORS !! %%
Neighbors_cellID(4).ON_Par_neighbor_id  = [302 631 7681] ; 
Neighbors_cellID(4).control_id    = []; %[ 2027 3736];
%}