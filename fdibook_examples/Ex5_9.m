% Example 5.9 - Solution of an AFDP using AFDSYN

% Uses the Control System Toolbox, Descriptor Systems Tools (DSTOOLS V0.61) 
% and Fault Detection and Isolation Tools (FDITOOLS V0.85)

clear variables
% define s as an improper transfer function
s = tf('s');
% define Gu(s), Gw(s), Gf(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     % enter Gu(s)
Gw = [1/(s+2); 0];                   % enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           % enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       % enter dimensions

% build model with additive faults 
sysf = fdimodset(ss([Gu Gf Gw]),struct('c',1:mu,'f',mu+(1:mf),'n',mu+mf+(1:mw)));

% design Q as [Q1;Q2]: enforce result using only the first basis vector of 
% the initial left nullspace basis [ I -Gu ] to design Q1 
tol = 1.e-7;    % tolerance for rank computations
opt_afdsyn = struct('tol',tol,'minimal',true,'smarg',-3,'poles',[-2,-3],...
                    'rdim',2,'nullspace',false,'HDesign',[1 0]);
[Q,R,info] = afdsyn(sysf,opt_afdsyn);  

% permute rows of Q and R to match example
ip = [2 1];
Q = minreal(tf(Q(ip,:)),tol) 
Rf = minreal(tf(R(ip,'faults')),tol)
Rw = minreal(tf(R(ip,'noise')),tol)

% display the resulting gap
info.gap

% this gap corresponds to the minimum 
min(fdif2ngap(R,[],[info.S;info.S2]))