% Example 5.13 - Solution of an EMMP using Procedure EMMS
% Uses the Control System Toolbox and Descriptor System Tools

% define s as an improper transfer function
s = tf('s');
% enter Gu(s), Gf(s) and Mr(s)
Gu = [s/(s^2+3*s+2) 1/(s+2);
     s/(s+1) 0;
      0 1/(s+2)];
Gf = [s/(s^2+3*s+2) 1/(s+2);
     s/(s+1) 0;
      0 1/(s+2)];
Mr = tf(eye(2));                  % enter Mr(s)
[p,mf] = size(Gf); mu = size(Gu,2);

% compute left nullspace basis as Q1(s) = [ I -Gu(s) ]
Q1 = [eye(p) -Gu]; Rf = Gf;

% % alternative computation
% opt_gln = struct('tol',1.e-7,'m2',mf);
% Q_Rf = glnull(ss([Gu Gf; eye(mu,mu+mf)]),opt_gln); 
% Q1 = Q_Rf(:,1:p+mu); Rf = Q_Rf(:,p+mu+1:end);

% check solvability condition
if rank(evalfr(Rf,rand)) ~= mf
   error('No solution exist')
end

% check for unstable or infinite zeros
gzero(ss(Rf))   % zeros at infinity and in the origin exist 

tol = 1.e-7;                     % set tolerance
sdeg = -1;                       % set stability degree

% using minimal dynamic covers, compute Q2 such that Q2*Q1 is a lower order
% left annihilator and Q2*Rf full row rank; 
% select basis vectors [2 3] of Q1 to combine them with basis vector 1
% to obtain Rf_Q = Q2*[Rf Q1], with Q2*Rf in the leading mf columns and 
% Q2*Q1 in the trailing p+mu columns
cinp = [ 2 3 1];
Rf_Q = glmcover1(ss([Rf(cinp,:) Q1(cinp,:) ]),mf,tol);


% compute the irreducible realization of Qtilde = Mr*(inv(Rf)*Q) by 
% first solving the linear rational matrix equation Rf*X = Q
X = grsol(Rf_Q,p+mu,struct('tol',tol));
Qtilde = gir(Mr*X,tol);

% compute stable and proper Q = Q4*Qtilde with suitable diagonal Q4 = M
Q = ss(zeros(0,p+mu)); M = ss(zeros(0,0));
opt_glcf = struct('tol',tol,'sdeg',sdeg);
for i=1:mf
    [Qi,Mi] = glcf(Qtilde(i,:),opt_glcf);
    sc = get(zpk(Mi),'k');   % scale with gain to fit example
    Q = [Q;Qi/sc]; M = append(M,Mi/sc);
end
% convert to standard state space representation
Q = gss2ss(Q); M = gss2ss(M); 
% display results
minreal(tf(Q)), tf(M)

% check solution
G = [Gu Gf;eye(mu,mu+mf)]; F = [zeros(mf,mu) M*Mr];
norm(Q*G-F,inf)


