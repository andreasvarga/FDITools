% Example 5.7 - Solution of an AFDP using AFDSYN

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

% computation using [ I -Gu ] as initial left nullspace basis 
% enforce using the sum of basis vectors via a design matrix H1 = [1 1]
tol = 1.e-7;    % tolerance for rank computations
opt_afdsyn = struct('tol',tol,'minimal',false,'smarg',-2,'sdeg',-2,...
    'rdim',1,'nullspace',false,'HDesign',[1 1],'nonstd',2,'epsreg',0.2);
[Q,R,info] = afdsyn(sysf,opt_afdsyn); info 
% Q = -Q; R = -R;   % scale with -1 to match example
Q = minreal(tf(Q),tol) 
Rf = minreal(tf(R(:,'faults')),tol)
Rw = minreal(tf(R(:,'noise')),tol)
gap = info.gap

