% Example 5.16 - Solution of an H-infinity AMMP
% Uses the Control System Toolbox and Descriptor System Tools

% define system with control, noise and actuator fault inputs 
A = [-.8 0 0; 0 -.5 .6; 0 -.6 -.5];
Bu = [1 1;1 0;0 1]; Bw = 0.25*[0 0;0 1;1 0]; Bf = Bu;
C = [0 1 1; 1 1 0]; Du = zeros(2,2);
% define Gu(s), Gw(s), Gf(s) and Mr(s)
Gu = ss(A,Bu,C,0); Gw = ss(A,Bw,C,0); Gf = Gu;
mf = size(Gf,2); Mr = ss(eye(mf)); 
[p,mu] = size(Gu); mw = size(Gw,2); n = order(Gu); 

% compute left nullspace basis as Q1(s) = [ I -Gu(s) ]
Q1 = ss(A,[zeros(n,p) -Bu],C,[eye(p) -Du]); Rf = Gf; Rw = Gw;

% check solvability condition
if rank(evalfr(Rf,rand)) ~= mf
   error('No solution')
end

% check for unstable or infinite zeros of [Rf Rw]
Rf_Rw = ss(A,[Bf Bw],C,0);
gzero(Rf_Rw)                 % two infinite zeros

tol = 1.e-7;                 % set tolerance
sdeg = -10;                  % set stability degree

% compute the quasi-co-outer-co-inner factorization of [Rf Rw]
[Gi,Go] = goifac(Rf_Rw,struct('tol',tol)); 

% compute Q = inv(Go)*Q1 using explicit formulas 
Qbar = dss([Go.a Go.b; Go.c Go.d],[Q1.b; Q1.d],...
        [ zeros(mf,n) -eye(mf)], zeros(mf,p+mu),...
        [eye(n,n+mf); zeros(mf,n+mf)]);

% compute [F1 F2 ] = [Mr 0]*Gi'
F1_F2 = [Mr zeros(mf,mw)]*Gi'; 

% solve L-inf LDP
options = struct('tol',tol,'reltol',5.e-4); 
[Q4,mininf] = glinfldp(F1_F2,mw,options);
% update Q to make it stable and proper
Qtilde = Q4*Qbar; 
Q = ss(zeros(0,p+mu)); M = ss(zeros(0,0));
opt_glcf = struct('tol',tol,'sdeg',sdeg,'mindeg',true,'mininf',true);
for i=1:mf
    [Qi,Mi] = glcf(Qtilde(i,:),opt_glcf);
    % normalize Mi to unit H-infinity norm to match example
    sc = norm(Mi,inf)*sign(dcgain(Mi));  
    Q = [Q;Qi/sc]; M = append(M,Mi/sc);
end

% convert solution to standard state space
Q = gss2ss(gir(Q,tol)); M = gss2ss(M); 

% check solution
G = [Gf Gw Gu; zeros(mu,mf+mw) eye(mu)]; 
F = M*Mr*eye(mf,mu+mw+mf);
err_sub = norm(Q*G-F,inf)

% compare with the (improper) optimal solution
Yopt = glinfldp(M*F1_F2,mw,options);
Qopt = Yopt*Qbar;
err_opt = norm(gir(Qopt*G-F,tol),inf)

%% simulate parametric step responses 
figure

% setup system with additive actuator faults 
sysftest = ss(A,[Bu Bf],C,0);
% set input groups
set(sysftest,'InputGroup',struct('controls',1:mu,'faults',mu+(1:mf)));

% form internal form of the filter R := Q*[Gu Gf; I 0];
R = Q*[sysftest;eye(mu,mu+mf)]; 

% visualization of step responses of all inputs 
% set names for residual and input components
R.OutputName = strcat(strseq('r_{',1:size(Q,1)),'}');
R.InputName = [strcat(strseq('u_{',1:mu),'}');
               strcat(strseq('f_{',1:mf),'}')];

% simulate step responses from all inputs
step(R(:,{'faults','controls'}),10) 
drawnow
hold on

% choose a ngrid-by-ngrid grid for parameter values
ngrid = 5;
rhovals = -0.25:.5/(ngrid-1):0.25;
for i=1:ngrid 
    for j=1:ngrid
       rho1 = rhovals(i); rho2 = rhovals(j);
       Ap = [-.8 0 0;0 -.5*(1+rho1) .6*(1+rho2); 0 -0.6*(1+rho2) -0.5*(1+rho1)];
       sysftest.a = Ap;
       R = Q*[sysftest;eye(mu,mu+mf)]; 
       step(R(:,{'faults','controls'}),10) 
       drawnow
    end
end
title('Step Responses')
ylabel('Residuals')

