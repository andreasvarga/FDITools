% Example 5.5 - Solution of an AFDP using EFDSYN

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

% use Q = [ I -Gu ] as initial left nullspace basis 
tol = 1.e-7;    % tolerance for rank computations
opt_efdsyn = struct('tol',tol,'minimal',false,'smarg',-2,'sdeg',-3,'nullspace',false);
[Q,R,info] = efdsyn(sysf,opt_efdsyn); info 
Q = minreal(tf(Q),tol) 
Rf = minreal(tf(R(:,'faults')),tol)
Rw = minreal(tf(R(:,'noise')),tol)
gap = fdif2ngap(R)
% the gap can be arbitraly increased by setting the stability margin and stabiity degree! 

