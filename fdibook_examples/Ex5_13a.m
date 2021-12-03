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

tol = 1.e-7;                     % set tolerance
sdeg = -1;                       % set stability degree

% G = [Gu Gf;eye(mu,mu+mf)]; F = [zeros(mf,mu) Mr];
G = [Gu Gf;eye(mu,mu+mf)]; F = [zeros(mf,mu) Mr];
Qtilde = glsol(ss(G),ss(F),struct('tol',tol,'mindeg',true)); 
norm(gminreal(Qtilde*G-F,tol),inf)

% compute stable and proper Q = Q4*Qtilde with suitable diagonal Q4 = M
Q = ss(zeros(0,p+mu)); M = ss(zeros(0,0));
opt_glcf = struct('tol',tol,'sdeg',sdeg);
for i=1:mf
    [Qi,Mi] = glcf(Qtilde(i,:),opt_glcf);
    sc = get(zpk(Mi),'k');   % scale to have the same M as in example
    Q = [Q;Qi/sc]; M = append(M,Mi/sc);
end
% convert to standard state space representation
Q = gss2ss(Q); M = gss2ss(M); minreal(tf(Q)), tf(M)

% check solution
G = [Gu Gf;eye(mu,mu+mf)]; F = [zeros(mf,mu) M*Mr];
norm(gss2ss(gir(Q*G-F,tol),tol),inf)

% Note: For the same Mr, the solution Q is different from that 
% computed in Example 5.13!


