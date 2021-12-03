% Example 5.10 - Solution of an EFDIP using EFDISYN
% Uses the Control System Toolbox, DSTOOLS V0.6 (and later), and 
% the FDITOOLS V0.2 (and later)

% enter output and fault vector dimensions
p = 3; mf = 3;
rng('default')
nu = floor(1+4*rand); mu = floor(1+4*rand);
nd = floor(1+4*rand); md = floor(1+4*rand);
% define random Gu(s) and Gd(s) with triplex sensor redundancy
% and Gf(s) = I for triplex sensor faults
Gu = ones(3,1)*rss(nu,1,mu); % enter Gu(s) in state-space form
Gd = ones(3,1)*rss(nd,1,md); % enter Gd(s) in state-space form

% build synthesis model with sensor faults
sysf = fdimodset([Gu Gd],struct('c',1:mu,'d',mu+(1:md),'fs',1:3));

SFDI = [ 0 1 1; 1 0 1; 1 1 0] > 0; % enter structure matrix

% set options for least order synthesis with EFDISYN
options = struct('tol',1.e-7,'sdeg',-1,'rdim',1,'SFDI',SFDI); 
[Qt,Rft] = efdisyn( sysf, options );

% normalize Q and Rf to match example
scale = sign([ Rft{1}.d(1,2) Rft{2}.d(1,3) Rft{3}.d(1,1)]);
for i = 1:3, Qt{i} = scale(i)*Qt{i}; Rft{i} = scale(i)*Rft{i}; end
Q = vertcat(Qt{:}); Rf = vertcat(Rft{:});

% check synthesis conditions: Q[Gu Gd;I 0] = 0 and Q[Gf; 0] = Rf
syse = [sysf;eye(mu,mu+md+mf)]; % form Ge = [Gu Gd Gf;I 0 0];
norm_Ru_Rd = norm(Q*syse(:,{'controls','disturbances'}),inf)
norm_rez = norm(Q*syse(:,'faults')-Rf,inf)

% check weak and strong fault detectability
S_weak = fditspec(Rf)
[S_strong,abs_dcgains] = fdisspec(Rf)

% evaluate step responses
set(Rf,'InputName',strseq('f_',1:mf),'OutputName',strseq('r_',1:size(SFDI,1)));
step(Rf);
title('Step responses from the fault inputs'), ylabel('')
