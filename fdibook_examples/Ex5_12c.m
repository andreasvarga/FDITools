% Example 5.12 - Solution of an EMMP using EMMSYN
% Uses the Control System Toolbox, DSTOOLS V0.71 and 
% the FDITOOLS V1.0 

% enter output and fault vector dimensions
p = 3; mf = 3;
% generate random dimensions for system order and input vectors
rng('default')
nu = floor(1+4*rand); mu = floor(1+4*rand);
nd = floor(1+4*rand); md = floor(1+4*rand);
% define random Gu(s) and Gd(s) with triplex sensor redundancy
% and Gf(s) = I for triplex sensor faults
Gu = ones(3,1)*rss(nu,1,mu); % enter Gu(s) in state-space form
Gd = ones(3,1)*rss(nd,1,md); % enter Gd(s) in state-space form

% build synthesis model with sensor faults
sysf = fdimodset([Gu Gd],struct('c',1:mu,'d',mu+(1:md),'fs',1:3));

% enter reference model for the TFM from faults to residual
Mr = fdimodset(ss([ 0 1 -1; -1 0 1; 1 -1 0]),struct('f',1:mf));

% solve an exact model-matching problem using EMMSYN
[Q,R,info] = emmsyn(sysf,Mr);

% check the synthesis: Q*Ge = M*Me and R = M*Mr, where
% Ge = [Gu Gd Gf; I 0 0] and Me = [0 0 Mr ].
Ge = [sysf; eye(mu,mu+md+mf)]; Me = [zeros(p,mu+md) Mr];
norm(gminreal(Q*Ge)-info.M*Me,inf)
norm(R-info.M*Mr,inf)

% compare with the one step solution
% solve Qbar*Ge = Me 
Qbar = glsol(Ge,Me,struct('tol',1.e-7));

% compare solutions by computing the infinity norm of Q-Qbar
norm(Q-Qbar,inf)