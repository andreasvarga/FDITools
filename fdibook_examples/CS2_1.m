%% CS2_1  - Case-study example: Monitoring air data sensor faults
%           Robust least order LTI synthesis
clear variables
%% Part 1 - Model setup
% load aircraft multiple-model SYSACSM and 
% percentages of mass variations massi
load('cs2data.mat')
% build minimal realizations of AC-models with sensor faults
% set dimensions
%SYSACSM = SYSACSM([1:4 7,8],:);
[p,m,N] = size(SYSACSM); nom = 7;             % nominal system 
% set sensor indices and set dimensions of inputs
%     [ AoA VCAS ]
sen = [  2    4  ]; 
md = 2; mu = m-md; mf = length(sen); n = max(order(SYSACSM));

% form system with faults [Gu Gd Gf] 
idm = eye(p);
syssenf = [ SYSACSM idm(:,sen)];  % add fault inputs
% set input groups
syssenf.InputGroup.controls = 1:mu;          % controls
syssenf.InputGroup.disturbances = mu+(1:md); % disturbances
syssenf.InputGroup.faults = mu+md+(1:mf);    % faults

Mr = eye(mf);                       % set desired reference model

%% Part 2 - Setup of the synthesis specifications

% compute the achievable strong specifications for constant faults
opt = struct('tol',1.e-7,'FDTol',1.e-5,'FDGainTol',0.001,...
    'FDFreq',0,'sdeg',-0.05,'m1',mu+md);
% apply genspec to [Gu Gd Gf; I 0 0]
S_strong = fdigenspec([syssenf(:,:,nom); eye(mu,mu+md+mf)],opt)

% define SFDI, the signatures for isolation of single faults
SFDI = eye(mf) > 0;


%% Part 3 - Multiple filter synthesis using Procedure EFDI
% set options for least order synthesis with EFDISYN
options = struct('tol',1.e-7,'sdeg',-5,'smarg',-0.05,'simple',false,...
    'FDFreq',0,'FDGainTol',0.0001,'rdim',1,'minimal',true,'SFDI',SFDI);
clear Q, Q(:,:,1:N) = ss(zeros(mf,p+mu));
for i = 1:N
   Qi = efdisyn( syssenf(:,:,i), options );  % synthesis of Qi
   Q(:,:,i) = [Qi{1};Qi{2}];
end

% % alternative computation 
% Q(:,:,1:N) = ss(zeros(mf,p+mu));
% for i = 1:N, 
%     cde = [syse(:,:,i).c syse(:,:,i).d]; 
%     Q(:,:,i) = ss([zeros(2,9) eye(2)]/cde); 
% end
% % 

%%  Part 4 - Assesment of synthesis results for individual syntheses
% form extended system [ Gu Gd Gf; I 0 0] with outputs extended with u
syse = [syssenf;eye(mu,mu+md+mf)];     

% form [Ru Rd Rf] 
R = Q*syse;

% check robustness by computing the H-inf norms in grid points
Nref = norm(R-[zeros(mf,mu+md) eye(2)],inf)

% visualization of step responses of all inputs 
% set names for residual and input components
set(R,'OutputName',{'r_1','r_2'}) 
for i=1:mu, R.InputName{i}=['u_{',num2str(i),'}']; end
for i=1:md, R.InputName{mu+i}=['d_{',num2str(i),'}']; end
for i=1:mf, R.InputName{mu+md+i}=['f_{',num2str(i),'}']; end
figure(1), step(R,10), grid, % simulate step responses from all inputs
ylabel(''), title('')

%%  Part 5 - Assesment of synthesis results for nominal synthesis
% form extended system [ Gu Gd Gf; I 0 0] with outputs extended with u
syse = [syssenf;eye(mu,mu+md+mf)];     

% form [Ru Rd Rf] 
R = Q(:,:,nom)*syse; 

% check robustness by computing the H-inf norms in grid points
Nnom = norm(R-[zeros(mf,mu+md) eye(2)],inf)

% visualization of step responses of all inputs 
% set names for residual and input components
set(R,'OutputName',{'r_1','r_2'}) 
for i=1:mu, R.InputName{i}=['u_{',num2str(i),'}']; end
for i=1:md, R.InputName{mu+i}=['d_{',num2str(i),'}']; end
for i=1:mf, R.InputName{mu+md+i}=['f_{',num2str(i),'}']; end
figure(2), step(R,10), grid, % simulate step responses from all inputs
ylabel(''), title('')

%%  Part 6 - Assesment of synthesis results for constant interpolated gain

% perform gain fitting using 0th order interpolation of elements
Q0 = zeros(mf,p+mu);
for i = 1:mf
    for j = 1:p+mu
        Q0(i,j) = polyfit(massi(:),squeeze(Q(i,j,:).d),0);
    end
end

% form extended system [ Gu Gd Gf; I 0 0] with outputs extended with u
syse = [syssenf;eye(mu,mu+md+mf)];     

% form [Ru Rd Rf] 
R = Q0*syse; 

% check robustness by computing the H-inf norms in grid points
N0 = norm(R-[zeros(mf,mu+md) eye(2)],inf)

% visualization of step responses of all inputs 
% set names for residual and input components
set(R,'OutputName',{'r_1','r_2'}) 
for i=1:mu, R.InputName{i}=['u_{',num2str(i),'}']; end
for i=1:md, R.InputName{mu+i}=['d_{',num2str(i),'}']; end
for i=1:mf, R.InputName{mu+md+i}=['f_{',num2str(i),'}']; end
figure(3), step(R,10), grid, % simulate step responses from all inputs
ylabel(''), title('')

%% Part 7 - Multiple model synthesis of a constant gain 

% define parameterized constant filter gain 
Qbar = ltiblock.ss('Qbar',0,mf,p+mu); 
Qbar.d.Value = Q(:,:,nom).d; 
Qbar.d.Free(1,4) = false;  % enforce VCAS decoupling 
Qbar.d.Free(2,2) = false;  % enforce AoA decoupling 

% define soft objective
syse = [syssenf;eye(mu,mu+md+mf)];          
E = (Qbar*syse-[zeros(mf,mu+md) eye(mf)]); 
E.InputName = 'udf'; E.OutputName = 'r';
Soft = TuningGoal.Gain('udf','r',1); 

% perform optimal tuning
Eopt = systune(E,Soft,[]);
Qbaropt = getValue(Qbar, Eopt); 

% scale to unit DC-gains
sc=dcgain(Qbaropt*syse(:,mu+md+1:end,nom));
Qbaropt = sc\Qbaropt 

% form [Ru Rd Rf] 
R = Qbaropt*syse; 

% check robustness by computing the H-inf norms in grid points
Nopt = norm(R-[zeros(mf,mu+md) eye(2)],inf)

% visualization of step responses of all inputs 
% set names for residual and input components
set(R,'OutputName',{'r_1','r_2'}) 
for i=1:mu, R.InputName{i}=['u_{',num2str(i),'}']; end
for i=1:md, R.InputName{mu+i}=['d_{',num2str(i),'}']; end
for i=1:mf, R.InputName{mu+md+i}=['f_{',num2str(i),'}']; end
figure(4), step(R,10), grid, % simulate step responses from all inputs
ylabel(''), title('')


[Nref Nnom N0 Nopt]


