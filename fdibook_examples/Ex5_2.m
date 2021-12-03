% Example 5.2 - EFDP has no solution 
% Uses the Control System Toolbox, DSTOOLS V0.61 (and later), and 
% the FDITOOLS V0.85 

% define s as an improper transfer function
s = tf('s');
% define Gu(s), Gd(s), Gf(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     % enter Gu(s)
Gd = [1/(s+2); 0];                   % enter Gd(s)
Gf = [(s+1)/(s+2) 0; 0 1];           % enter Gf(s)
p = 2; mu = 1; md = 1; mf = 2;       % set dimensions

% build model with additive faults 
sysf = fdimodset(ss([Gu Gd Gf]),struct('c',1:mu,'d',mu+(1:md),'f',mu+md+(1:mf)));

% call of EFDSYN with default options 
try
  [Q,Rf] = efdsyn(sysf); 
catch err
  disp(err.message)
  % compute achievable structure matrix
  S_achievable = fdigenspec(sysf)  
end
 