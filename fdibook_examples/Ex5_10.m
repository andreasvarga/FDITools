% Example 5.10 - Solution of an EFDIP
% Uses the Control System Toolbox and Descriptor System Tools

% enter output and fault vector dimensions
p = 3; mf = 3;
% generate random dimensions for system order and input vectors
nu = floor(1+4*rand); mu = floor(1+4*rand);
nd = floor(1+4*rand); md = floor(1+4*rand);
% define random Gu(s) and Gd(s) with triplex sensor redundancy
% and Gf(s) for triplex sensor faults
Gu = ones(3,1)*rss(nu,1,mu); % enter Gu(s) in state-space form
Gd = ones(3,1)*rss(nd,1,md); % enter Gd(s) in state-space form
Gf = eye(3);                 % enter Gf(s) for sensor faults
tol = 1.e-7;                 % tolerance for rank tests

% build model with faults
sysf = [Gu Gd Gf];         
% set input groups
sysf.InputGroup.controls = 1:mu;          % controls
sysf.InputGroup.disturbances = mu+(1:md); % disturbances
sysf.InputGroup.faults = mu+md+(1:mf);    % faults

S = [ 0 1 1; 1 0 1; 1 1 0];  % enter structure matrix

%% Step 1) of Procedure EFDI
% compute Q1, the left nullspace of [Gu Gd;I 0], and  Rf1 = Q1*[Gf;0]
% QRf contains Q1 in QRf(:,1:p+mu) and Rf1 in QRf(:,p+mu+1:p+mu+mf)
options_glnull=struct('tol',tol,'m2',mf); 
QRf = glnull([sysf; eye(mu,mu+md+mf)],options_glnull); 

%% Step 2): of Procedure EFDI 
% initialize overall filter Q and fault detection system Rf
nb = size(S,1);      % number of filters $n_b$
Qt = cell(nb,1); Rft = cell(nb,1);
% options for EFDSYN for the synthesis of scalar output filters
options = struct('tol',tol,'rdim',1); 
QRf.InputGroup.aux = 1:p+mu+mf ;         % indices of [Q1 Rf1]
for i = 1:nb
    % Step 2.1): define  Rf11 = Rf1(:,indd), Rf12 = Rf1(:,indf)
    indd = find(S(i,:) == 0); indf = find(S(i,:) ~= 0);
    QRf.InputGroup.disturbances = p+mu+indd; % indices of Rf11
    QRf.InputGroup.faults = p+mu+indf;       % indices of Rf12
    
    % Step 2.2): apply  Procedure EFD to {0,Rf11,Rf12,[Q1 Rf1]}
    % to determine a least order Q1i such that Q1i*Rf11 = 0
    % QRfauxi contains: [Q1i*Rf12 Q1i*Q1 Q1i*Rf1]
    [~,QRfauxi] = efdsyn( QRf, options ); 
    QRfi = QRfauxi(:,'aux');  % extract [Q1i*Q1 Q1i*Rf1]
    QRfi.InputGroup.aux = []; % clear auxiliary input group 
    
    % Step 2.3): extract Q{i} = Q1i*Q1 and Rf{i} = Q1i*Rf1
    Qt{i} = QRfi(:,1:p+mu);
    Rft{i} = QRfi(:,p+mu+(1:mf));
end
%%

% normalize Q and Rf to match example
scale = sign([ Rft{1}.d(1,2) Rft{2}.d(1,3) Rft{3}.d(1,1)]);
for i = 1:3, Qt{i} = scale(i)*Qt{i}; Rft{i} = scale(i)*Rft{i}; end
Q = [Qt{1};Qt{2};Qt{3}], Rf = [Rft{1};Rft{2};Rft{3}] 
