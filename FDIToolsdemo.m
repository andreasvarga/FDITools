% FDITOOLSDEMO    Demonstration of Fault Detection and Isolation Filter 
%                 Synthesis Tools (FDITOOLS).
%                 A. Varga.

%  Copyright 2015-2021 A. Varga
%
clc
format compact
echo on
% FDITOOLS is a collection of functions for the analysis and solution 
% of fault detection problems. The FDITOOLS collection relies on the 
% the Control System Toolbox and the Descriptor System Tools (DSTOOLS). 
%
% The functions implemented in FDITOOLS are based on the computational 
% procedures for the synthesis of fault detection filters described in 
% the book (V,2017):
%
%   Varga A. 
%   Solving Fault Diagnosis Problems - Linear Synthesis Techniques, 
%   volume 84 of Studies in Systems, Decision and Control, 
%   Springer International Publishing, xxii+382, 2017. 

pause % Press any key to continue ...    

% The current release of FDITOOLS is version 1.0.6, dated May 1, 2021.
% The available funtions in FDITOOLS are:
%
% Setup of synthesis models.
%   fdimodset  - Setup of models for solving FDI synthesis problems.
%   mdmodset   - Setup of models for solving model detection synthesis problems.
%
% FDI related analysis.
%   fdigenspec - Generation of achievable FDI specifications.
%   fdichkspec - Feasibility analysis of a set of FDI specifications.
%
% Model detection related analysis.
%   mddist     - Computation of distances between component models. 
%   mddist2c   - Computation of distances to a set of component models. 
%
% Performance evaluation of FDI filters
%   fditspec   - Computation of the weak or strong structure matrix.
%   fdisspec   - Computation of the strong structure matrix.
%   fdifscond  - Fault sensitivity condition of FDI filters. 
%   fdif2ngap  - Fault-to-noise gap of FDI filters.
%   fdimmperf  - Model-matching performance of FDI filters.
%
% Performance evaluation of model detection filters
%   mdperf     - Distance mapping performance of model detection filters.
%   mdmatch    - Distance matching performance of model detection filters.
%   mdgap      - Noise gaps of model detection filters.
%   
% Synthesis of FDI filters.
%   efdsyn     - Exact synthesis of fault detection filters.
%   afdsyn     - Approximate synthesis of fault detection filters.
%   efdisyn    - Exact synthesis of fault detection and isolation filters.
%   afdisyn    - Approximate synthesis of fault detection and isolation filters.
%   emmsyn     - Exact model matching based synthesis of FDI filters.
%   ammsyn     - Approximate model matching based synthesis of FDI filters.
%
% Synthesis of model detection filters.
%   emdsyn     - Exact synthesis of model detection filters.
%   amdsyn     - Approximate synthesis of model detection filters.

pause % Press any key to continue ...  

%%
clc
% Example 1 - Solution of an exact fault detection problem (EFDP) 
clear variables

% We use the unstable model in Example 5.4 of (V,2017) 
s = tf('s');                         % define complex variable s
% define Gu(s), Gd(s), Gf(s)
Gu = [(s+1)/(s-2); (s+2)/(s-3)];     % enter Gu(s)
Gd = [(s-1)/(s+2); 0];               % enter Gd(s)
p = 2; mu = 1; md = 1;  mf = 2;      % set dimensions

% setup the synthesis model with faults
sysf = fdimodset(ss([Gu Gd]),struct('c',1,'d',2,'f',1,'fs',2));

pause % Press any key to continue ...  

clc
% call of EFDSYN with the options for stability degree -3 and the synthesis 
% of a scalar output filter
[Q,Rf] = efdsyn(sysf,struct('sdeg',-3,'rdim',1)); 

% display the implementation form Q and internal form Rf of the 
% resulting fault detection filter
tf(Q), tf(Rf)

pause % Press any key to continue ...  

clc
% check synthesis results: R(s) := Q(s)*Ge(s) = [ 0  0  Rf(s)], 
% with Ge(s) = [Gu(s) Gd(s)  Gf(s);I  0 0]
norm_rez = fdimmperf(Q * [sysf;eye(mu,mu+md+mf)],Rf)

pause % Press any key to continue ...  

% check weak and strong fault detectability
S_weak = fditspec(Rf)
[S_strong,abs_dcgains] = fdisspec(Rf)

pause % Press any key to continue ...  

% determine the fault sensitivity condition
FSCOND = fdifscond(Rf,0)

pause % Press any key to continue ...  

% check synthesis using step responses from the faults
Rf.OutputName = 'r'; Rf.InputName = strseq('f_',1:size(Rf,2));
step(Rf)
title('Step responses from the fault inputs'), ylabel('')

pause % Press any key to continue ...  


%% 
clc
% Example 2 - Solution of an exact fault detection and isolation problem (EFDIP) 
clear variables

% We build a random system with triplex sensor faults
% as in  Example 5.12 of (V,2017) 

% set the output and fault vector dimensions
p = 3; mf = 3;
% generate random dimensions for system orders, control and disturbance inputs
rng('default')
nu = floor(1+4*rand); mu = floor(1+4*rand);
nd = floor(1+4*rand); md = floor(1+4*rand);
% define random Gu(s) and Gd(s) with triplex sensor redundancy
% and Gf(s) = I for triplex sensor faults
Gu = ones(3,1)*rss(nu,1,mu); % enter Gu(s) in state-space form
Gd = ones(3,1)*rss(nd,1,md); % enter Gd(s) in state-space form
Gf = eye(3);                 % enter Gf(s) for sensor faults

% build synthesis model with sensor faults
sysf = fdimodset([Gu Gd],struct('c',1:mu,'d',mu+(1:md),'fs',1:3));

pause % Press any key to continue ...  

% determine the achievable strong specifications for constant faults
opt = struct('tol',1.e-7,'FDTol',1.e-5,'FDGainTol',0.01,...
             'FDFreq',0,'sdeg',-0.05);
% apply fdigenspec to [Gu Gd Gf; I 0 0]
S_strong = fdigenspec(sysf,opt)

pause % Press any key to continue ...  
clc
% choose specifications for a voting based FDI system
SFDI = S_strong(sum(S_strong,2)==2,:)

% set options for least order synthesis with EFDISYN
options = struct('tol',1.e-7,'rdim',1,'SFDI',SFDI,'minimal',true); 

% call of EFDISYN to determine a bank of three scalar output filters
[Qt,Rft] = efdisyn( sysf, options );

pause % Press any key to continue ...  

% normalize Q(s) and Rf(s) to match example
sc = sign([Rft{1}.d(1,2) Rft{2}.d(1,3) Rft{3}.d(1,1)]);
for i = 1:3
    Qt{i} = sc(i)*Qt{i}; Rft{i} = sc(i)*Rft{i}; 
end
Q = vertcat(Qt{:}); Rf = vertcat(Rft{:});
inpnames = [strseq('y',1:p');strseq('u',1:mu)]; 
outnames = strseq('r',1:size(SFDI,1)); 
fnames = strseq('f',1:mf);
Q  = set(Q,'InputName',inpnames,'OutputName',outnames)
Rf = set(Rf,'InputName',fnames,'OutputName',outnames)           

pause % Press any key to continue ...  

clc
% check synthesis results: R(s) := Q(s)*Ge(s) = [ 0  0  Rf(s)],
% with Ge(s) = [Gu(s) Gd(s)  Gf(s);I  0 0]
R = gir(Q * [sysf;eye(mu,mu+md+mf)]);
norm_rez = fdimmperf(R,Rf)

pause % Press any key to continue ...  

% check achieved strong structure matrix
% check fault isolability
isequal(SFDI,fdisspec(Rf))

pause % Press any key to continue ...  

% determine achieved fault sensitivity conditions
FSCOND = fdifscond(Rf,0,SFDI)

pause % Press any key to continue ...  

% check synthesis using step responses from the faults
figure
Rf.OutputName = strseq('r_',1:size(SFDI,1));
Rf.InputName = strseq('f_',1:mf);
step(Rf,5), ylabel('')

pause % Press any key to continue ...  

%%
clc
% Example 3 - Solution of an exact model matching problem (EMMP) 

% We use the input data and the results of the previous example 
% to set up an equivalent EMMP (as in  Example 5.12 of (V,2017) )

% we use as reference model the resulting global filter Rft
Mr = Rft;
pause % Press any key to continue ...  

% we solve an EMMP with EMMSYN using the default option for 
% least order synthesis 
[Q,R,info] = emmsyn(sysf,Mr);
% the updating factor M{i} is contained in info.M{i}
tf(vertcat(info.M{:}))

pause % Press any key to continue ...  

clc
% check synthesis results: R(s) := Q(s)*Ge(s) = M(s)*[0 0 Mr(s)],
% with Ge(s) = [Gu(s) Gd(s) Gf(s); I 0 0]
Ge = [sysf;eye(mu,mu+md+mf)]; Rtemp = vertcat(Q{:})*Ge;
norm_rez = fdimmperf(Rtemp,append(info.M{:})*vertcat(Mr{:}))

pause % Press any key to continue ...  

% check computed R(s)
norm_Rdif = fdimmperf(Rft,R)  

pause % Press any key to continue ...  

% compare the resulting filter Q with the previously computed filter Qt
norm_Qdif = fdimmperf(Q,Qt,inf)
pause % Press any key to continue ...  

%%
clc
% Example 4 - Solution of a strong EFDIP 
clear variables

% We use Example 5.13 of (V,2017)
s = tf('s');                         % define complex variable s
% define Gu(s)
Gu = [s/(s^2+3*s+2) 1/(s+2);         % enter Gu(s)
     s/(s+1) 0;
      0 1/(s+2)];
[p,mu] = size(Gu); mf = mu;          % set dimensions

% build model with faults
sysf = fdimodset(ss(Gu),struct('c',1:mu,'f',1:mu));

% define Mr(s)
Mr = fdimodset(ss(eye(2)),struct('f',1:mf));

pause % Press any key to continue ...  

clc
% call of EMMSYN with the options for stability degree -1 and 
% tolerances for rank computations and observability tests set to 1.e-7 
emmsyn_options = struct('sdeg',-1,'tol',1.e-7,'tolmin',1.e-7);
[Q,R,info] = emmsyn(sysf,Mr,emmsyn_options); 
% the updating factor M is contained in info.M
tf(info.M)
pause % Press any key to continue ...  

clc
% verify the model matching condition: Q*Ge = M(s)*[0 Mr(s)],
% with Ge(s) = [Gu(s) Gf(s); I 0]
Ge = [sysf;eye(mu,mu+mf)]; Rtemp = Q*Ge;
norm_rez = fdimmperf(Rtemp,info.M*Mr)

pause % Press any key to continue ...  

% check computed R(s)
norm_dif = fdimmperf(Rtemp,R)  

pause % Press any key to continue ...  


% verify reference model updating: R = M*Me
norm_rdif = norm(R-info.M*Mr,inf)

pause % Press any key to continue ...  

%%
clc
% Example 5 - Solution of an exact model detection problem (EMDP) 
clear variables

% We use Example 6.1 of (V,2017), representing the 
% lateral dynamics model (without faults) of an F-16 aircraft with:
% n  = 4 states
% mu = 2 control inputs
% p  = 4 measurable outputs

A = [-.4492 0.046 .0053 -.9926;
       0    0     1     0.0067;
   -50.8436 0   -5.2184  .722;
    16.4148 0     .0026 -.6627];
Bu = [0.0004 0.0011; 0 0; -1.4161 .2621; -0.0633 -0.1205]; 
C = eye(4); p = size(C,1); mu = size(Bu,2); 

pause % Press any key to continue ...  

% define the loss of efficiency (LOE) faults as input scaling gains
% Gamma(i,:) = [ 1-rho1(i) 1-rho2(i) ]
Gamma = 1 - [ 0  0 0 .5 .5 .5 1  1 1;
              0 .5 1  0 .5  1 0 .5 1 ]';           
N = size(Gamma,1);  % number of LOE cases

pause % Press any key to continue ...  

% define a multiple physical fault model Gui = Gu*diag(Gamma(i,:))
sysu = ss(zeros(p,mu,N,1));
for i = 1:N
    sysu(:,:,i,1) = ss(A,Bu*diag(Gamma(i,:)),C,0);
end

% set input groups
set(sysu,'InputGroup',struct('controls',1:mu));

pause % Press any key to continue ...  
clc
% call of EMDSYN with the options for stability degree -1 and pole -1 for 
% the filters, tolerance and a design matrix H to form a linear combination
% of the left nullspace basis vectors
H = [ 0.7645 0.8848 0.5778 0.9026 ];
emdsyn_options = struct('sdeg',-1,'poles',-1,'HDesign',H);

pause % Press any key to continue ...  

[Q,R,info] = emdsyn(sysu,emdsyn_options); 

% inspect achieved performance
info.MDperf

pause % Press any key to continue ...  

clc
% plot the step responses for the internal filter representations
figure
echo off
k1 = 0;
for j = 1:N
  k1 = k1+1; 
  k = k1;
  for i=1:N 
    subplot(N,N,k), 
    [r,t] = step(R{j,i},4); 
    plot(t,r(:,:,1),t,r(:,:,2)), 
    if i == 1, title(['Model ',num2str(j)]), end
    if i == j, ylim([-1 1]), end 
    if j == 1, ylabel(['r^(^', num2str(i),'^)'],'FontWeight','bold'), end
    if i == N && j == 5, xlabel('Time (seconds)','FontWeight','bold'), end
    k = k+N;
  end
end

%%
echo on
pause % Press any key to continue ...  

clc
% Example 6 - Separation of disturbance inputs from noise inputs}
clear variables

% We use the model obtained by recasting parametric uncertainties as noise
% inputs in Example 2.2 of (V,2017)

% define system with control, noise and  fault inputs
A = [-.8 0 0; 0 -.5 .6; 0 -.6 -.5];
Bu = [1 1;1 0;0 1]; Bw = [0 0;0 1;1 0]; 
C = [0 1 1; 1 1 0];
Du = zeros(2,2); Dw = Du; 
mu = 2; mw = 2; mf = 2; % set input dimensions, two actuator faults

pause % Press any key to continue ...  

% setup the model with noise inputs redefined as faults
sysf = fdimodset(ss(A,[Bu Bw],C,[Du Dw]),struct('c',1:2,'f',1:4));

pause % Press any key to continue ...  

% compute the achievable weak specifications
S_weak = fdigenspec(sysf)
% zeros in the leading mw columns indicate possible exact decoupling of
% some of noise inputs

pause % Press any key to continue ...  

% select specifications containing nonzero elements in the 
% trailing mf columns
S_red = S_weak(sum(S_weak(:,end-mf+1:end),2) > 0,:)

pause % Press any key to continue ...  

% select specifications containing zeros in the leading mw elements
for i = 1:mw
    S_red(S_red(:,i) == 0,:)
end

pause % Press any key to continue ...  


%%
clc
% Example 7 - Solution of an approximate fault detection problem
clear variables

% We use Example 5.3 of (V,2017) (modified)
s = tf('s'); % define the Laplace variable s
Gu = [(s+1)/(s+2); (s+2)/(s+3)];  % enter Gu(s)
Gw = [(s-1)/(s+2); 0];            % enter Gw(s)

% build the model with additive faults having Gf(s) = [Gu(s) eye(p)];
sysf = fdimodset(ss([Gu Gw]),struct('c',1,'f',1,'fs',1:2,'n',2));
mu = 1; mw = 1; p = 2; mf = mu+p; % set dimensions

pause % Press any key to continue ...  

clc
% try first an exact synthesis based design using EFDSYN
[Q,R] = efdsyn(sysf);  % R = [ Rf  Rw ]

pause % Press any key to continue ...  

% compute the achieved initial gap
gap0 = fdif2ngap(R) 

pause % Press any key to continue ...  

clc
% perform synthesis with  AFDSYN
[Q,R,info] = afdsyn(sysf);  % R = [ Rf  Rw ]

tf(Q), tf(R(:,'faults')),   tf(R(:,'noise'))

% the gap has been improved 
[info.gap gap0] 

pause % Press any key to continue ...  

clc
% check results: R(s) := Q(s)*Ge(s) = [0 Rf(s) Rw(s)],
% with Ge(s) = [Gu(s) Gf(s) Gw(s); I 0 0]
norm_rez = fdimmperf(Q*[sysf;eye(mu,mu+mf+mw)],R)

pause % Press any key to continue ...  

% check strong fault detectability
S_strong = fdisspec(R)

pause % Press any key to continue ...  

% simulate step responses from fault and noise inputs
figure
inpnames = {'f_1','f_2','f_3','w'};
set(R,'InputName',inpnames,'OutputName','r');
step(R); ylabel('')
title('Step responses from fault and noise inputs')

pause % Press any key to continue ...  


%%
clc
% Example 8 - Solution of an approximate fault detection and isolation problem
clear variables

% We use Example 5.3 of (V,2017) (modified)
s = tf('s'); % define the Laplace variable s
Gu = [(s+1)/(s+2); (s+2)/(s+3)];  % enter Gu(s)
Gw = [(s-1)/(s+2); 0];            % enter Gw(s)

% build the model with additive faults having Gf(s) = [Gu(s) eye(p)];
sysf = fdimodset(ss([Gu Gw]),struct('c',1,'f',1,'fs',1:2,'n',2));
mu = 1; mw = 1; p = 2; mf = mu+p; % set dimensions

pause % Press any key to continue ...  


% select the desired structure matrix SFDI
S = fdigenspec(sysf); SFDI = S(sum(S,2)==2,:)
nb = size(SFDI,1);

pause % Press any key to continue ...  

clc
% try first an exact synthesis based design using EFDISYN

% set options for least order synthesis with EFDISYN
options = struct('tol',1.e-7,'smarg',-3,'sdeg',-3,'SFDI',SFDI);
[Q,R] = efdisyn(sysf,options);

pause % Press any key to continue ...  

% compute achieved initial gaps 
format short e
gap0 = fdif2ngap(R,[],SFDI) 

pause % Press any key to continue ...  

clc
% perform an approximate synthesis with AFDISYN

% set options for least order synthesis with AFDISYN
options = struct('tol',1.e-7,'smarg',-3,'sdeg',-3,'SFDI',SFDI);
[Q,R,info] = afdisyn(sysf,options);  % R{i} = [ Rf{i} Rw{i} ]
pause % Press any key to continue ...  

format short e
% gaps have been (partly) improved
[gap0 info.gap] 

pause % Press any key to continue ...  

clc
% check R{i}(s) := Q{i}(s)*Ge(s) = [0 Rf{i}(s) Rw{i}(s)],
% with Ge(s) = [Gu(s) Gd(s) Gf(s); I 0 0]
Ge = [sysf;eye(mu,mu+mf+mw)]; 
Rtemp = cell(nb,1);
for i = 1:nb, Rtemp{i} = Q{i}*Ge; end
norm_rez = fdimmperf(Rtemp,R)

pause % Press any key to continue ...  

% check fault isolability of constant faults
S_strong = fdisspec(R) 
match_SFDI = isequal(SFDI,S_strong)

pause % Press any key to continue ...  

clc
figure 
% simulate step responses from fault and noise inputs
inpnames = {'f_1','f_2','f_3','w'};
outnames = {'r_1','r_2','r_3'};
Rtot = [R{1};R{2};R{3}];
set(Rtot,'InputName',inpnames,'OutputName',outnames);
step(Rtot); ylabel('Residuals')
title('Step responses from fault and noise inputs')

pause % Press any key to continue ...  

%%
clc
% Example 9 - Solution of an approximate model-matching problem
clear variables

% This is Example 5.11 of (V,2017)
s = tf('s'); % define the Laplace variable s
Gu = [(s+1)/(s+2); (s+2)/(s+3)];  % enter G_u(s)
Gf = [(s+1)/(s+2) 0; 0 1];        % enter G_f(s)
Gw = [1/(s+2); 0];                % enter G_w(s)
mu = 1; mf = mu+1; mw = 1; p = 2; % set dimensions

% build the synthesis model with additive faults
inputs = struct('c',1:mu,'f',mu+(1:mf),'n',mu+mf+(1:mw));
sysf = fdimodset(ss([Gu Gf Gw]),inputs);
pause % Press any key to continue ...  

% determine maximal achievable structure matrix
tol = 1.e-7; % set tolerance for rank computations
Smax = fdigenspec(sysf,struct('tol',tol))

% choose the targeted reference model
Mr = fdimodset(ss(eye(mf)),struct('faults',1:mf));
pause % Press any key to continue ...  

clc
% set options for AMMSYN
options = struct('tol',tol,'reltol',1.e-8);

% perform an approximate model-matching based synthesis with AMMSYN
[Q,R,info] = ammsyn(sysf,Mr,options);

% optimal and suboptimal performance and achieved gaps
info
format short e
gap = fdif2ngap(R,[],fditspec(Mr))

pause % Press any key to continue ...  

clc
% check decoupling condition and synthesis results
Ge = [sysf;eye(mu,mu+mf+mw)]; Rtemp = gminreal(Q*Ge);
norm_rez  = fdimmperf(Rtemp,R)

pause % Press any key to continue ...  

% check achieved synthesis performance 
gammaopt = fdimmperf(Rtemp,Mr)

pause % Press any key to continue ...  

%%
clc
% Example 10 - Solution of an approximate model detection problem (AMDP) 
clear variables

% We use the model of Example 6.2 of (V,2017), representing the 
% lateral dynamics model of an F-16 aircraft with:
% n  = 4 states
% mu = 2 control inputs
% p  = 2 measurable outputs
% mw = 6 noise inputs 
A = [-.4492 .046  .0053 -.9926;
       0    0     1      .0067;
   -50.8436 0   -5.2184  .722;
    16.4148 0     .0026 -.6627];
Bu = [.0004 .0011; 0 0; -1.4161 .2621; -.0633 -.1205];
[n,mu] = size(Bu); p = 2; mw = n+p; m = mu+mw;
Bw = eye(n,mw);
C = 180/pi*eye(p,n);  Du = zeros(p,mu); Dw = [zeros(p,n) eye(p)]; 

pause % Press any key to continue ...  

% define the loss of efficiency (LOE) faults as input scaling gains
% Gamma(i,:) = [ 1-rho1(i) 1-rho2(i) ]
Gamma = 1 - [ 0  0 0 .5 .5 .5 1  1 1;
              0 .5 1  0 .5  1 0 .5 1 ]';           
N = size(Gamma,1);  % number of LOE cases

pause % Press any key to continue ...  

% define a multiple physical fault model with the i-th model 
% [ Gui Gw ], where Gui = Gu*diag(Gamma(i,:))
sysuw = ss(zeros(p,mu+mw,N,1));
for i = 1:N
    sysuw(:,:,i,1) = ss(A,[Bu*diag(Gamma(i,:)) Bw],C,[Du Dw]); 
end

% setup synthesis model
sysuw = mdmodset(sysuw,struct('controls',1:mu,'noise',mu+(1:mw)));

pause % Press any key to continue ...  

clc
% use preliminary (nonminimal) design with  EMDSYN
opt_emdsyn = struct('sdeg',-1,'poles',-1,'minimal',false);
[Q,R,info0] = emdsyn(sysuw,opt_emdsyn); 

% achieved performance
info0.MDperf

pause % Press any key to continue ...  

clc
% compute achieved initial gaps
format short e
gap0 = mdgap(R)

pause % Press any key to continue ...  

clc
% use nonminimal design with  AMDSYN
% set options for AMDSYN
opt_amdsyn = struct('sdeg',-1,'poles',-1,'minimal',false);

% perform an approximate synthesis with AMDSYN
[Q,R,info] = amdsyn(sysuw,opt_amdsyn); info.MDperf

% achieved initial and optimal gaps 
[gap0, info.MDgap]
%

pause % Press any key to continue ...  

% plot performance measures
figure
mesh(info.MDperf)
colormap hsv
title('Norms of final residual models')
ylabel('Residual numbers')
xlabel('Model numbers')

pause % Press any key to continue ...  


% generate input signals for Ex. 6.2
% set signal amplitude scaling
d = diag([ 1 1 0.01 0.01 0.01 0.01 0.03 0.03]);
t = (0:0.01:10)';  ns = length(t);
% generate sinusoidal and square-wave inputs
usin = gensig('sin',pi,t(end),0.01)+1.5;
usquare = gensig('square',pi*2,t(end),0.01)+0.3;
% use zero mean white noise inputs 
u = [ usquare usin (rand(ns,mw)-0.5)]*d;  
 
%
clc
% plot the step responses for the internal filter representations
figure
echo off
k1 = 0; alpha = 0.9; beta = 0.1; gamma = 10;
for j = 1:N,
  k1 = k1+1; k = k1;
  for i=1:N, 
    subplot(N,N,k), 
    [r,t] = lsim(R{j,i},u,t);
    % use a Narendra filter with (alpha,beta,gamma) = (0.9,0.1,10)
    theta = alpha*sqrt(r(:,1).^2+r(:,2).^2)+...
        beta*sqrt(lsim(tf(1,[1 gamma]),r(:,1).^2+r(:,2).^2,t));
    plot(t,real(theta)), 
    if i == 1, title(['Model ',num2str(j)]), end
    if i == j, ylim([0 4]), end 
    if j == 1, ylabel(['\theta_', num2str(i)],'FontWeight','bold'), end
    if i == N && j == 5, xlabel('Time (seconds)','FontWeight','bold'), end
    k = k+N;
  end
end

pause % Press any key to continue ...  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format
