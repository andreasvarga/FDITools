%% CS1_1  - Case-study example: Monitoring flight actuator faults
%           No measurements of control surface angles are used
clear variables
%% Part 1 - Model setup
% load aircraft multiple-model SYSACM, actuator model ACT,
% output-feedback gain K, percentages of mass variations massi
load('cs1data.mat')
% build minimal realizations of AC-models with actuator faults
% set dimensions
[p,m,N] = size(SYSACM); nom = 7;             % nominal system 
% set primary actuator indices
%           [ ailerons,  elevators, stabilizer, ruder ]
act_prim = [ [1,2,15,16], [17,19],     18,       20   ]; 
mu = size(ACT,1); md = m-mu; mf = length(act_prim);

% form system with faults [Gu Gd Gf] and set dimensions
sysact = SYSACM*append(ACT,eye(md)); % connect actuators
sysact = minreal(sysact);            % build minimal realizations
sysactf = sysact(:,[1:m act_prim]);  % add fault inputs
% set input groups
sysactf.InputGroup.controls = 1:mu;          % controls
sysactf.InputGroup.disturbances = mu+(1:md); % disturbances
sysactf.InputGroup.faults = mu+md+(1:mf);    % faults

% determine closed-loop stability margin
sdeg_cloop = max(max(real(eig(feedback(sysact,K,1:mu,1:p,+1)))))

%% Part 2 - Setup of the synthesis specifications
% compute the achievable weak specifications
opt = struct('tol',1.e-7,'FDTol',1.e-5,'m1',mu+md);
% apply genspec to [Gu Gd Gf; I 0 0]
S = fdigenspec([sysactf(:,:,nom); eye(mu,mu+md+mf)],opt);
%
% compute the achievable strong specifications for constant faults
opt = struct('tol',1.e-7,'FDTol',1.e-5,'FDGainTol',0.01,...
    'FDFreq',0,'sdeg',-0.05,'m1',mu+md);
% apply genspec to [Gu Gd Gf; I 0 0]
S_strong = fdigenspec([sysactf(:,:,nom); eye(mu,mu+md+mf)],opt)


% define SFDI, the signatures for isolation of single faults
SFDI = [0 1 1 1 1 0 1 0
        1 0 1 1 0 1 1 0
        1 1 0 1 1 0 1 0
        1 1 1 0 1 1 0 0
        1 1 1 1 0 0 0 0
        0 0 0 0 0 0 0 1] > 0;

%% Part 3 - Synthesis using Procedure EFDI 
% set options for least order synthesis with EFDISYN
options = struct('tol',1.e-7,'sdeg',-5,'smarg',-0.5,'simple',false,...
    'FDFreq',0,'FDGainTol',0.0001,'rdim',1,'SFDI',SFDI);
[Q,Rf,info] = efdisyn( sysactf(:,:,nom), options );  

%%  Part 4 - Assesment of synthesis results
% form extended system [ Gu Gd Gf; I 0 0] with outputs extended with u
syse = [sysactf;eye(mu,mu+md+mf)];     
% form extended closed-loop system with output feedback u = K*y+v
sysefb = feedback(syse,K,1:mu,1:p,+1); 

%% build overall detector and fault detection system for open-loop system
Qtot = ss(zeros(0,p+mu)); Rftilde = ss(zeros(0,mf));
for i = 1:size(SFDI,1)
    Qtot = [Qtot; Q{i}];
    Rftilde = [Rftilde; Rf{i}];
end

% open-loop checks: Q[Gu Gd;I 0] = 0$ and Q[Gf; 0] = Rftilde
norm(gir(Qtot*syse(:,1:mu+md,nom),5.e-5),inf)
norm(gir(Qtot*syse(:,mu+md+1:end,nom)-Rftilde,1.e-5),inf)

% check of achieved structure matrix
if size(Rftilde,1) == size(SFDI,1)
   if any(any((abs(dcgain(Rftilde)) > .5) - SFDI))
      error(['Desired FDI specification is not feasible'])
   end
end

% evaluate [Ru Rd Rf] for the closed-loop setup
Rtot = Qtot*sysefb;
% check robustness by computing the H-inf norms 
NormRu = squeeze(norm(Rtot(:,'controls'),inf));
NormRd = squeeze(norm(Rtot(:,'disturbances'),inf));
NormRfmRfnom = squeeze(norm(Rtot(:,'faults')-Rftilde,inf));

[NormRu, NormRd, NormRfmRfnom]
%%
figure, plot(massi,NormRu,massi,NormRd,massi,NormRfmRfnom)

% visualization of step responses of fault inputs 
rdim = size(Rftilde,1); resn = cell(1,rdim); faultn = cell(1,mf); 
for i=1:rdim, resn{i} = ['r_{',num2str(i),'}']; end
set(Rftilde,'OutputName',resn)  % set names for residual components
for i=1:mf, faultn{i} = ['f_{',num2str(i),'}']; end
set(Rftilde,'InputName',faultn) % set names for fault components
figure, step(Rftilde,10), grid, % simulate step responses for faults

% set names for residual and fault components
set(Rtot,'OutputName',resn)  
for i=1:mf, Rtot.InputName{mu+md+i}=faultn{i}; end
figure, step(Rtot(:,'faults'),10), grid, % simulate step responses for faults

