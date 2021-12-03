% Example 5.16 - Solution of an H-infinity AMMP using function AMMSYN

% Uses the Control System Toolbox, Descriptor Systems Tools (DSTOOLS V0.71) 
% and Fault Detection and Isolation Tools (FDITOOLS V1.0)

clear variables
% define system with control and noise inputs 
A = [-.8 0 0; 0 -.5 .6; 0 -.6 -.5];
Bu = [1 1;1 0;0 1]; Bw = 0.25*[0 0;0 1;1 0]; 
C = [0 1 1; 1 1 0]; Du = zeros(2,2); Dw = zeros(2,2);
[p,mu] = size(Du); mf = mu; mw = size(Dw,2); 

% setup synthesis model with additive actuator faults with Bf = Bu, Df = Du
sysf = fdimodset(ss(A,[Bu Bw],C,[Du Dw]),struct('c',1:mu,'n',mu+(1:mw),'f',1:mu)); 

% define Mr(s) = I
Mr = fdimodset(ss(eye(2)),struct('f',1:mf));

options = struct('tol',1.e-7,'nullspace',false,'reltol',5.e-4,'sdeg',-10);
[Q,R,info] = ammsyn(sysf,Mr,options);
gamma_opt0 = info.gammaopt0  % optimal performance for the initial problem 
gamma_opt  = info.gammaopt   % optimal performance 
gamma_sub  = info.gammasub   % suboptimal performance 

% check suboptimal solution
Ge = [sysf;eye(mu,mu+mf+mw)]; Me = [zeros(mf,mu) Mr zeros(mf,mw)];
norm(R-Q*Ge(:,{'faults','noise'}),inf)

%% simulate parametric step responses 
figure

% setup system with additive actuator faults 
sysftest = fdimodset(ss(A,Bu,C,0),struct('c',1:mu,'f',1:mu));

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
