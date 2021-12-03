% Example 5.4 - Solution of an EFDP
% Uses the Control System Toolbox and Descriptor System Tools

% define s as an improper transfer function
s = tf('s');
% define Gu(s), Gd(s), Gf(s)
Gu = [(s+1)/(s-2); (s+2)/(s-3)];     % enter Gu(s)
Gd = [(s-1)/(s+2); 0];               % enter Gd(s)
Gf = [(s+1)/(s-2) 0; (s+2)/(s-3) 1]; % enter Gf(s)
p = 2; mu = 1; md = 1; mf = 2;       % set dimensions

% compute a left nullspace basis Q of [Gu Gd; I 0]
Q1 = glnull(ss([Gu Gd;eye(mu,mu+md)]));

% compute Rf1 = Q1[Gf;0]
Rf1 = gir(Q1*[Gf;zeros(mu,mf)]);

% check solvability using a random frequency
if min(abs(evalfr(Rf1,rand))) > 0.01
   % compute a stable left coprime factorization [Q1 Rf1]=inv(Q3)*[Q,Rf]
   % enforce stability degree -3
   [Q_Rf,Q3] = glcf([Q1,Rf1],struct('sdeg',-3));
   % extract Q and Rf
   Q = Q_Rf(:,1:p+mu); Rf = Q_Rf(:,p+mu+1:end); 
   % normalize Q and Rf to match example
   scale = evalfr(Rf(1),inf);
   Q = tf(Q/scale), Rf = tf(Rf/scale)
else
   disp('No solution exists')
end

