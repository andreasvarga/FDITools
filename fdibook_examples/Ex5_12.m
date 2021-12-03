% Example 5.12 - Solution of an EMMP using Procedure EMM
% Uses the Control System Toolbox and Descriptor System Tools

% enter output and fault vector dimensions
p = 3; mf = 3;
% generate random dimensions for system order and input vectors
nu = floor(1+4*rand); mu = floor(1+4*rand);
nd = floor(1+4*rand); md = floor(1+4*rand);
% define random Gu(s) and Gd(s) with triplex sensor redundancy
% and Gf(s) for triplex sensor faults
Gu = ones(3,1)*rss(nu,1,mu); % enter Gu(s) in state-space form
Gd = ones(3,1)*rss(nd,1,md); % enter Gd(s) in state-space form
Gf = eye(3);                 % enter Gf(s) for sensor faults
% enter reference model
Mr = ss([ 0 1 -1; -1 0 1; 1 -1 0]);

% two step solution using Procedure EMM
% 1. Compute left nullspace and the reduced system
Q_Rf = glnull([Gu Gd Gf; eye(mu,mu+md+mf)],struct('m2',mf));
Q1 =  Q_Rf(:,1:p+mu); Rf = Q_Rf(:,p+mu+1:end);

% 2. Solve X*Rf = Mr and update Q
Q2 = glsol(Rf,Mr,struct('tol',1.e-7));
Q = Q2*Q1;

% one step solution
% solve Qbar*Ge = Me, where Ge = [Gu Gd Gf; I 0 0] and Me = [0 0 Mr ].
Ge = [Gu Gd Gf; eye(mu,mu+md+mf)]; Me = [zeros(p,mu+md) Mr];
Qbar = glsol(Ge,Me,struct('tol',1.e-7));

% compare solutions by computing the infinity norm of Q-Qbar
norm(Q-Qbar,inf)