% Example 5.13 - Solution of a strong EMMP using EMMSYN
% Uses the Control System Toolbox, DSTOOLS V0.71  and 
% the FDITOOLS V1.0 

clear variables
% define s as an improper transfer function
s = tf('s');
% enter Gu(s)
Gu = [s/(s^2+3*s+2) 1/(s+2);
     s/(s+1) 0;
      0 1/(s+2)];
[p,mu] = size(Gu); mf = mu;

% build model with faults
sysf = fdimodset(ss(Gu),struct('c',1:mu,'f',1:mu));

% define Mr(s) = I
Mr = fdimodset(ss(eye(2)),struct('f',1:mf));

% solve a strong EFDIP using EMMSYN (for an invertible reference model) 
opts_emmsyn = struct('tol',1.e-7,'sdeg',-1,'minimal',false);
[Q,R,info] = emmsyn(sysf,Mr,opts_emmsyn);

% check solution
G = [sysf;eye(mu,mu+mf)]; 
norm(Q*G-info.M*[zeros(mf,mu) Mr],inf)

% display results
minreal(tf(Q)), tf(info.M)

