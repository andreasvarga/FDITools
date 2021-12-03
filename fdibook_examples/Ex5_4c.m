% Example 5.4 - Solution of an EFDP using EFDSYN
% Uses the Control System Toolbox, DSTOOLS V0.61 and 
% the FDITOOLS V0.85 

% define s as an improper transfer function
s = tf('s');
% define Gu(s) and Gd(s)
Gu = [(s+1)/(s-2); (s+2)/(s-3)];     % enter Gu(s)
Gd = [(s-1)/(s+2); 0];               % enter Gd(s)
Gf = [(s+1)/(s-2) 0; (s+2)/(s-3) 1]; % enter Gf(s)
p = 2; mu = 1; md = 1; mf = 2;       % set dimensions

% setup the synthesis model with faults with Gf(s) = [ Gu(s) [0;1]]
sysf = fdimodset(ss([Gu Gd]),struct('c',1,'d',2,'f',1,'fs',2));

% call of EFDSYN with the options for stability degree -3 and the synthesis 
% of a scalar output filter
[Q,Rf] = efdsyn(sysf,struct('sdeg',-3,'rdim',1)); 

% normalize Q and Rf to match example
scale = evalfr(Rf(1,1),inf);
Q = tf(Q/scale), Rf = tf(Rf/scale) 
