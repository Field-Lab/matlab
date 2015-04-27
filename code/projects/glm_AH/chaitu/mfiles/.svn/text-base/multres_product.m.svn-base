% Chaitu Ekanadham, 04/07/2009

% Compute the product AWA'x = b where 
% A = [Ac;Af] is blocked at 2 different resolutions
% W  = diag(wf) ""
% x = [xc;xf] ""

function b = multres_product(Ac,Af,wf,xc,xf)

r = size(Af,2)/size(Ac,2);

Acxc = Ac'*xc; %eval(sprintf('whos AC xc')) %coarse Tc x b % check
Afxf = Af'*xf; % fine Tf x b

WAcxc = mult_row(reprows(Acxc,r),wf); % fine (because of W) - Tf x b
WAfxf = mult_row(Afxf,wf); % fine Tf x b

% Block the fine one to coarse resolution

WAcxc_blocked = block_spikes(WAcxc',r)'; % Tc x b
WAfxf_blocked = block_spikes(WAfxf',r)'; % Tc x b

AcWAcxc = Ac*WAcxc_blocked; %eval(sprintf('whos WAcxc_blocked')) % check
AfWAcxc = Af*WAcxc;

AcWAfxf = Ac*WAfxf_blocked; %eval(sprintf('whos Ac WAfxf_blocked')) % check
AfWAfxf = Af*WAfxf;

b = [(AcWAcxc + AcWAfxf);(AfWAcxc + AfWAfxf)];

%fprintf('multres_product\n')
%keyboard
