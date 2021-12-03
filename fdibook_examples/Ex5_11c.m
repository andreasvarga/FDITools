% Example 5.11 - Solution of an AFDIP using AFDISYN

% Uses the Control System Toolbox, Descriptor Systems Tools (DSTOOLS V0.71) 
% and Fault Detection and Isolation Tools (FDITOOLS V1.0)

% define s as an improper transfer function
s = tf('s');
% define Gu(s), Gw(s), Gf(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     % enter Gu(s)
Gw = [1/(s+2); 0];                   % enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           % enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       % enter dimensions

% build the synthesis model with additive faults 
sysf = fdimodset(ss([Gu Gf Gw]),struct('c',1:mu,'f',mu+(1:mf),'n',mu+mf+(1:mw)));

% enter structure matrix
SFDI = eye(mf) > 0;                         

options = struct('tol',1.e-7,'smarg',-3,'sdeg',-3,'SFDI',SFDI);
[Q,R,info] = afdisyn(sysf,options);
% scale with -1 to match example
Q{2} = -Q{2}; R{2} = -R{2};

Q1 = minreal(tf(Q{1})), Rf1 = minreal(tf(R{1}(:,'faults'))), Rw1 = minreal(tf(R{1}(:,'noise')))  
Q2 = minreal(tf(Q{2})), Rf2 = minreal(tf(R{2}(:,'faults'))), Rw2 = minreal(tf(R{2}(:,'noise')))  
info.gap

% computation of achieved gaps
format short e
gaps = fdif2ngap(R,[],SFDI)
