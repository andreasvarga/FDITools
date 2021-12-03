%% CS2_2  - Case-study example: Monitoring air data sensor faults
%           Robust least order LPV synthesis
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
mu = m-2; md = 2; mf = length(sen); n = max(order(SYSACSM));

% form system with faults [Gu Gd Gf] 
idm = eye(p);
syssenf = [ SYSACSM idm(:,sen)];  % add fault inputs
% set input groups
syssenf.InputGroup.controls = 1:mu;          % controls
syssenf.InputGroup.disturbances = mu+(1:md); % disturbances
syssenf.InputGroup.faults = mu+md+(1:mf);    % faults

syse = [syssenf;eye(mu,mu+md+mf)];  % form extened system
Mr = eye(mf);                       % set desired reference model

%% Part 2 - Setup of the synthesis specifications

% compute the achievable strong specifications for constant faults
opt = struct('tol',1.e-7,'FDTol',1.e-5,'FDGainTol',0.001,...
    'FDFreq',0,'sdeg',-0.05,'m1',mu+md);
% apply genspec to [Gu Gd Gf; I 0 0]
S_strong = fdigenspec([syssenf(:,:,nom); eye(mu,mu+md+mf)],opt)

% define SFDI, the signatures for isolation of single faults
SFDI = eye(2) > 0;

%% Part 3 - LPV synthesis using 4th order polynomial gain interpolation

% form extended system
syse = [syssenf;eye(mu,mu+md+mf)];  

% determine constant filter gains 
q = zeros(mf,p+mu,N);  
for i = 1:N, 
    cde = [syse(:,:,i).c syse(:,:,i).d]; 
    q(:,:,i) = [zeros(mf,n+mu+md) eye(mf)]/cde; 
end

k = 4;     % set order of polynomial approximation
% determine element-wise the $k$th order polynomial approximations 
% D1+x*D2+x^2*D3+x^3*D4+x^4*D5
theta = zeros(mf,p+mu,k+1);
x = 2*(massi-0.5);  %  normalize to [-1,1]
for i = 1:mf
    for j = 1:p+mu
        theta(i,j,k+1:-1:1) = polyfit(x(:),squeeze(q(i,j,:)),k);
    end
end

% evaluate gains in gridpoints using Horner's rule
% D1+x*(D2+x*(D3+x*(D4+x*D5)))
DQk = zeros(mf,p+mu,N);
for i = 1:N
    Qval = theta(:,:,k+1);
    for j = k:-1:1
        Qval = x(i)*Qval+theta(:,:,j);
    end
    DQk(:,:,i) = Qval;
end

% form [Ru Rd Rf] 
R = DQk*syse; 

% check robustness by computing the H-inf norms in grid points
Nk = norm(R-[zeros(mf,mu+md) eye(2)],inf)

% visualization of step responses of all inputs 
% set names for residual and input components
set(R,'OutputName',{'r_1','r_2'}) 
for i=1:mu, R.InputName{i}=['u_{',num2str(i),'}']; end
for i=1:md, R.InputName{mu+i}=['d_{',num2str(i),'}']; end
for i=1:mf, R.InputName{mu+md+i}=['f_{',num2str(i),'}']; end
figure(6), step(R,10), grid, % simulate step responses from all inputs
ylabel(''), title('')

%% Part 4 - LPV synthesis using 3rd order polynomial gain interpolation
% setup of the LPV synthesis to use SYSTUNE

% parametrize LPV detector
domain = struct('mass',massi(:));
shapefcn = @(x) [x,x^2,x^3];  
% define gain as D0+x*D1+x^2*D2+x^3*D3  and  theta = [D0 D1 D2 D3]
Qd = tunableSurface('Qd',q(:,:,nom),domain,shapefcn);
Qd.Coefficients.free(1,2:p+mu:end) = false;
Qd.Coefficients.free(2,4:p+mu:end) = false;
Qd.Coefficients.free(1,4:p+mu:end) = false;
Qd.Coefficients.free(2,2:p+mu:end) = false;

% define soft objective
syse = [syssenf;eye(mu,mu+md+mf)];          
E0 = (Qd*syse-[zeros(mf,mu+md) eye(mf)]); 
E0.InputName = 'udf'; E0.OutputName = 'r';
Soft = TuningGoal.Gain('udf','r',1); 

% perform optimal tuning
E0opt = systune(E0,Soft,[]);
Qlpv = setData(Qd,E0opt.Blocks.Qd.Value)

% form [Ru Rd Rf] 
R = Qlpv*syse; 

% check robustness by computing the H-inf norms in grid points
Nlpvopt = norm(R-[zeros(mf,mu+md) eye(2)],inf)

% visualization of step responses of all inputs 
% set names for residual and input components
set(R,'OutputName',{'r_1','r_2'}) 
for i=1:mu, R.InputName{i}=['u_{',num2str(i),'}']; end
for i=1:md, R.InputName{mu+i}=['d_{',num2str(i),'}']; end
for i=1:mf, R.InputName{mu+md+i}=['f_{',num2str(i),'}']; end
figure(7), step(R,10), grid, % simulate step responses from all inputs
ylabel(''), title('')

[Nk Nlpvopt]

