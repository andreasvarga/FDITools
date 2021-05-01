function [Q,R,info] = efdisyn(sysf,options)
%EFDISYN Exact synthesis of fault detection and isolation filters
%  [Q,R,INFO] = EFDISYN(SYSF,OPTIONS) solves the exact fault detection and 
%  isolation problem (EFDIP) for a LTI system SYSF with additive faults and 
%  for a given structure matrix SFDI (specified via the OPTIONS structure).
%  Two banks of stable and proper filters Q{1}, ..., Q{nb} and 
%  R{1}, ..., R{nb} are computed, where nb is the number of rows of the 
%  structure matrix SFDI. The filters Q{1}, ..., Q{nb} are the 
%  solution of the EFDIP and R{1}, ..., R{nb) are the corresponding 
%  internal forms. The two banks of nb filters are stored in the 1 x nb 
%  cell arrays Q and R, respectively. 
%
%  The continuous- or discrete-time system SYSF must be given in a standard
%  or descriptor state-space form, which corresponds to the input-output 
%  form  
%
%       y = Gu*u + Gd*d + Gf*f + Gw*w ,
%
%  with the Laplace- or Z-transformed plant outputs y, control inputs u, 
%  disturbance inputs d, fault inputs f, and noise inputs w, and with 
%  Gu, Gd, Gf, and Gw the corresponding transfer-function matrices. 
%  The inputs u, d, f, and w of SYSF correspond to four input groups 
%  named, respectively, {'controls','disturbances','faults','noise'}.
%
%  The i-th filter Q{i}, determined in a standard state-space form, 
%  generates the i-th residual signal r_i and corresponds to the 
%  input-output (implementation) form
%
%            r_i = Q{i}*[ y ] = Qy_i*y + Qu_i*u .
%                       [ u ]
%
%  The inputs y and u of each of the resulting filter Q{i} are grouped in 
%  two input groups {'outputs','controls'}, respectively, and the output 
%  group {'residuals'} is defined for the residuals r_i. 
%
%  The i-th filter R{i}, determined in a standard state-space form, is the 
%  internal form of Q{i}, generates the i-th residual signal r_i, and 
%  corresponds to the input-output form
%
%       r_i = Ru_i*u + Rd_i*d + Rf_i*f + Rw_i*w ,
%
%  where 
%
%    R{i} := [ Ru_i Rd_i Rf_i Rw_i} ] = Q{i}*[ Gu Gd Gf Gw ]. 
%                                            [ I  0  0  0  ]
%
%  The solution of the EFDIP ensures that Ru_i = 0, Rd_i = 0, and Rf_i 
%  has the structure matrix SFDI(i,:), the i-th row of the given
%  structure matrix SFDI (specified via the option structure OPTIONS). 
%  The inputs f and w of each of the resulting filter R{i} are grouped 
%  in two input groups {'faults','noise'}, respectively, and 
%  the output group {'residuals'} is defined for the residuals r_i. 
%
%  The resulting filters Q{i} and R{i} have observable state-space 
%  realizations (AQi,BQi,CQi,DQi) and (AQi,BRi,CQi,DRi), respectively, and 
%  thus share the observable pairs (AQi,CQi). 
%
%  The OPTIONS structure allows to specify various user options, as
%  follows:
%  OPTIONS.SFDI      - desired structure matrix to solve the EFDIP
%                      (Default: [ 1 ... 1 ], i.e., solve an exact fault
%                                detection problem)
%  OPTIONS.tol       - relative tolerance for rank computations
%                      (Default: internally determined value)
%  OPTIONS.tolmin    - absolute tolerance for observability tests
%                      (Default: internally determined value)
%  OPTIONS.FDTol     - threshold for fault detectability checks
%                      (Default: 0.0001)
%  OPTIONS.FDGainTol - threshold for strong fault detectability checks
%                      (Default: 0.01)
%  OPTIONS.rdim      - a vector q, whose i-th component q(i) specifies the 
%                      desired number of residual outputs for Q{i} and R{i},
%                      or a scalar q, which specifies the same desired    
%                      number of residual outputs for all Q{i} and R{i} 
%                      (Default: [], in which case: 
%                                if OPTIONS.HDesign{i} is empty, then
%                                   q(i) = 1, if OPTIONS.minimal = true, or
%                                   q(i) is the dimension of the left  
%                                   nullspace which underlies the synthesis 
%                                   of Q{i} and R{i}, if 
%                                   OPTIONS.minimal = false;
%                                if OPTIONS.HDesign{i} is non-empty, then 
%                                   q(i) is the row dimension of the i-th  
%                                   design matrix contained in 
%                                   OPTIONS.HDesign{i})
%  OPTIONS.FDFreq    - vector of frequency values for strong detectability 
%                      checks (Default: [])
%  OPTIONS.smarg     - prescribed stability margin for the poles of the 
%                      filters Q{i} and R{i}
%                      (Default: -sqrt(eps) for a continuous-time system 
%                               and 1-sqrt(eps) for a discrete-time system) 
%  OPTIONS.sdeg      - prescribed stability degree for the poles of the 
%                      filters Q{i} and R{i}
%                      (Default:  -0.05 for a continuous-time system and
%                                  0.95  for a discrete-time system) 
%  OPTIONS.poles     - specifies a complex conjugated set of desired poles 
%                      within the stability domain to be assigned for the 
%                      filters Q{i} and R{i} (Default: []) 
%  OPTIONS.nullspace - specifies the proper nullspace basis option for the 
%                      initial step
%                      true  - use minimal proper basis (default); 
%                      false - use full-order observer based basis; 
%                              this option can be only used for a proper 
%                              system without disturbance inputs 
%  OPTIONS.simple    - specifies the option to employ simple proper 
%                      left nullspace bases for synthesis 
%                      true  - use simple bases; the orders of the  
%                              basis vectors used to determine Q{i} 
%                              are provided in INFO.degs{i}
%                      false - no simple basis computed (default) 
%  OPTIONS.minimal   - specifies the option to perform least order 
%                      synthesis of the filters Q{i} and R{i}
%                      true  - perform least order synthesis (default)
%                      false - no least order synthesis performed  
%  OPTIONS.tcond     - maximum alowed condition number of the employed 
%                      non-orthogonal transformations (Default: 1.e4).
%  OPTIONS.FDSelect  - index set of desired filters to be designed
%                      (Default: 1:nb)
%  OPTIONS.HDesign   - nb-dimensional cell array; OPTIONS.HDesign{i} is a 
%                      full row rank design matrix employed for the 
%                      synthesis of the i-th fault detection filter 
%                      (Default: [])
%
%  INFO is a structure containing additional information, as follows: 
%  INFO.tcond        - nb-dimensional vector; INFO.tcond(i) contains
%                      the maximum of the condition numbers of the employed 
%                      non-orthogonal transformation matrices to determine
%                      the i-th filter component Q{i}; a warning is 
%                      issued if any INFO.tcond(i) >= OPTIONS.tcond.
%  INFO.degs         - nb-dimensional cell array; INFO.degs{i} contains 
%                      the increasingly ordered degrees of a left minimal   
%                      polynomial nullspace basis, which can serve 
%                      for the synthesis of the i-th filter component Q{i} 
%  INFO.HDesign      - nb-dimensional cell array; INFO.HDesign{i} is the 
%                      design matrix H_i employed for the synthesis of 
%                      the i-th fault detection filter.
%
%  See also EFDSYN. 

%  Copyright 2015-2018 A. Varga
%  Author:    A. Varga, 14-12-2015.
%  Revisions: A. Varga, 15-02-2016, 28-05-2016, 17-11-2016, 24-08-2017,
%                       06-07-2018.
%
%  Method: The synthesis Procedure EFDI from [1] is implemented, which 
%  relies on the nullspace-based synthesis method proposed in [2]. For more 
%  details on the least order synthesis of fault detection filters see [3].
%
%  References:
%  [1] Varga A.:
%      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%      Springer Verlag, 2017; sec. 5.4.
%  [2] Varga, A.:
%      On designing least order residual generators for fault detection 
%      and isolation. 
%      16th International Conference on Control Systems and  
%      Computer Science, Bucharest, Romania, 2007.
%  [3] Varga, A.:
%      On computing least order fault detectors using rational 
%      nullspace bases. 
%      IFAC SAFEPROCESS'03 Symposium, Washington DC, USA, 2003.

narginchk(1,2)
nargoutchk(0,3)

% check input system form
if ~isa(sysf,'ss')
   error('The input system SYSF must be an SS object')
end

% decode input information
if isfield(sysf.InputGroup,'controls')
    % controls
    inpu = sysf.InputGroup.controls;
    mu = length(inpu);  
else
    inpu = []; mu = 0;
end
if isfield(sysf.InputGroup,'disturbances')
    % disturbances
    inpd = sysf.InputGroup.disturbances;
    md = length(inpd);  
else
    inpd = []; md = 0;
end
if isfield(sysf.InputGroup,'faults')
    % faults
    inpf = sysf.InputGroup.faults;
    mf = length(inpf);  
else
    inpf = []; mf = 0;
end
if isfield(sysf.InputGroup,'noise')
    % noise
    inpw = sysf.InputGroup.noise;
    mw = length(inpw);  
else
    inpw = []; mw = 0;
end
 
m = mu+md+mf+mw;
p = size(sysf,1);

if nargin < 2
    options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

discr = (sysf.Ts > 0);  % system type (continuous- or discrete-time)
if discr
    smax = 1-sqrt(eps); smin = 0;
else
    smax = -sqrt(eps);  smin = -inf;
end
 
% decode options

% FDI structure matrix
if isfield(options,'SFDI')
   SFDI = options.SFDI;   
   validateattributes(SFDI, {'logical'},{'binary'},'','OPTIONS.SFDI') 
   if size(SFDI,2) ~= mf
      error('Structure matrix incompatible with the number of faults')
   end
else
   SFDI = true(1,mf);   % solve by default an EFDP
end
nb = size(SFDI,1);      % number of filters 

% tolerance for rank tests
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end

% tolerance for observability tests
if isfield(options,'tolmin')
   tolmin = options.tolmin;
   validateattributes(tolmin, {'double'},{'real','scalar','>=',0},'','OPTIONS.tolmin') 
else
   tolmin = 0;
end

% threshold for fault detectability checks
if isfield(options,'FDTol')
   FDTol = options.FDTol;
   validateattributes(FDTol, {'double'},{'real','scalar','>=',0},'','OPTIONS.FDTol') 
else
   FDTol = 0.0001;
end

% threshold for strong fault detectability checks
if isfield(options,'FDGainTol')
   FDGainTol = options.FDGainTol;
   validateattributes(FDGainTol, {'double'},{'real','scalar','>=',0},'','OPTIONS.FDGainTol') 
else
   FDGainTol = 0.01;
end

% desired number of filter outputs
if isfield(options,'rdim')
   rdim = options.rdim;
   if ~isempty(rdim)
      validateattributes(rdim, {'double'},{'integer','vector','>',0},'','OPTIONS.rdim')
      if length(rdim) == 1
         rdim = rdim*ones(1,nb);
      else
         if nb ~= length(rdim)
            error('The dimension of the vector OPTIONS.rdim must be equal to the row dimension of OPTIONS.SFDI')
         end
      end
   end
else
   rdim = [];
end

% frequency values for strong detectability checks
if isfield(options,'FDFreq')
   FDFreq = options.FDFreq;
   if ~isempty(FDFreq)
      validateattributes(FDFreq, {'double'},{'real','vector','>=',0},'','OPTIONS.FDFreq')
   end
else
   FDFreq = [];
end

% desired stability degree
if isfield(options,'sdeg')
   sdeg = options.sdeg;
   if ~isempty(sdeg)
      validateattributes(sdeg, {'double'},{'real','scalar','<=',smax,'>=',smin},'','OPTIONS.sdeg') 
   end
else
   sdeg = [];
end
if isempty(sdeg)
   % set desired stability degree
   if discr
      sdeg = 0.95;
   else
      sdeg = -0.05;
   end
end    

% stability margin
if isfield(options,'smarg')
   smarg = options.smarg;
   if ~isempty(smarg)
      validateattributes(smarg, {'double'},{'real','scalar','<=',smax,'>=',smin},'','OPTIONS.smarg') 
   end
else
   smarg = [];
end
if isempty(smarg)
   % set stability margin
   if discr
      smarg = 1-sqrt(eps);
   else
      smarg = -sqrt(eps);
   end
end    

% desired poles
if isfield(options,'poles')
   poles = options.poles;
   if ~isempty(poles)
      validateattributes(poles, {'double'},{'vector'},'','OPTIONS.poles')
   end
else
   poles = [];
end

% check that poles are stable and form a complex conjugate set
if ~isempty(poles)
    if ~isequal(sort(poles(imag(poles)>0)),sort(conj(poles(imag(poles)<0)))) 
        error('POLES must be a self-conjugated complex vector')
    end 
    if (discr && ~isempty(find(abs(poles) > 1-sqrt(eps), 1))) || ...
       (~discr && ~isempty(find(real(poles) > -sqrt(eps), 1))) 
          error('The elements of POLES must lie in the stability region')
    end
end    


% option for nullspace basis
if isfield(options,'nullspace')
   nullspace = options.nullspace; 
   validateattributes(nullspace, {'logical'},{'binary'},'','OPTIONS.nullspace') 
else
   nullspace = true;
end

% option for simple basis
if isfield(options,'simple')
   simple = options.simple; 
   validateattributes(simple, {'logical'},{'binary'},'','OPTIONS.simple') 
else
   simple = false;
end


% option for least order synthesis
if isfield(options,'minimal')
   minimal = options.minimal; 
   validateattributes(minimal, {'logical'},{'binary'},'','OPTIONS.minimal') 
else
   minimal = true;
end

% upper margin for condition  number of used transformations
if isfield(options,'tcond')
   tcond = options.tcond;
   validateattributes(tcond, {'double'},{'real','scalar','>=',1},'','OPTIONS.tcond') 
else
   tcond = 1.e4;
end

% selected filters
if isfield(options,'FDSelect')
   FDSelect = options.FDSelect;
   if isempty(FDSelect)
      FDSelect = 1:nb;
   else
      validateattributes(FDSelect, {'double'},{'vector','integer','>=',1,'<=',nb,'increasing'},'','OPTIONS.FDSelect')
   end
else
   FDSelect = 1:nb;
end

% imposed design option to form linear combinations of basis vectors
if isfield(options,'HDesign')
   HDesign = options.HDesign;
   if ~isempty(HDesign)
      validateattributes(HDesign, {'cell'},{'vector','numel',nb},'','OPTIONS.HDesign')
      for i = 1:nb
          if ~isempty(HDesign{i})
              if ~isempty(rdim) && size(HDesign{i},1) ~= rdim(i)
                 error(['Row dimension of OPTIONS.HDesign{',num2str(i),'} must be equal to OPTIONS.rdim{',num2str(i),'}'])
              end
              if size(HDesign{i},1) ~= rank(HDesign{i})
                 error(['OPTIONS.HDesign{',num2str(i),'} must have full row rank'])
              end
          end
      end
   end
else
   HDesign = cell(1,nb);
end

% Step 1) of Procedure EFDI
% compute Q1, the left nullspace of [Gu Gd;I 0], and  
% R = [Rf1 Rw1] = Q1*[Gf Gw;0 0]
% QR contains Q1 in QR(:,1:p+mu) and R in QR(:,p+mu+1:end),
% where Rf1 is in QR(:,p+mu+1:p+mu+mf)
% form [ Gu Gd Gf Gw; I 0 0 0 ] 

% set options for initial nullspace computation (use no simple basis)
opts_glnull = struct('tol',tol,'m2',mf+mw);
% syse = [sysf(:,[inpu inpd inpf inpw]); eye(mu,m)];
% QR = glnull(syse,opts_glnull); 

% Step 1): nullspace based reduction
% compute Q1, the left nullspace of [Gu Gd;I 0], and  
% R = [Rf1 Rw1] = Q1*[Gf Gw;0 0]
% QR contains Q1 in QR(:,1:p+mu) and R in QR(:,p+mu+1:end),
% where Rf1 is in QR(:,p+mu+1:p+mu+mf)
% form [ Gu Gd Gf Gw; I 0 0 0 ] 
%
if nullspace || md || rcond(sysf.e) < 1.e-7 
   % form [ Gu Gd Gf Gw; I 0 0 0] 
   syse = [sysf(:,[inpu inpd inpf inpw]); eye(mu,m)];
   %
   % compute a left nullspace basis Q = Q1 of G1 = [Gu Gd; I 0] = 0 and
   % obtain QR = [ Q R ], where R = [ Rf Rw ] = Q*[Gf Gw;0 0]
   QR = glnull(syse,opts_glnull); 
else
   % compute minimal basis as Q = Q1 = [ I -Gu] and set
   % QR = [ Q R ], where R = [ Gf Gw Ga ]
   QR = [ eye(p) dss(sysf.a,[-sysf.b(:,inpu) sysf.b(:,[inpf inpw])],...
          sysf.c,[-sysf.d(:,inpu) sysf.d(:,[inpf inpw])],sysf.e,sysf.Ts)]; 
end

% check solvability conditions
if size(QR,1) == 0,
   error('Empty nullspace basis: the EFDIP is not solvable')
end

% Step 2): of Procedure EFDI 
% initialize overall filters Q and R
Q = cell(nb,1); 
if nargout > 1
   R = cell(nb,1); 
end
if nargout > 2
   info.tcond = ones(nb,1); 
   info.degs = cell(nb,1);
   info.HDesign = cell(nb,1);
end
% options for EFDSYN for the synthesis of fault detection filters
optionsi = struct('tol',tol,'sdeg',sdeg,'smarg',smarg,'FDFreq',FDFreq, ...
    'FDTol',FDTol,'FDGainTol',FDGainTol,'rdim',[],'simple',simple, ...
    'minimal',minimal,'tolmin',tolmin,'tcond',tcond,'poles',poles,'HDesign',[]); 
QR.InputGroup.aux = 1:p+mu+mf+mw;       % input indices of [Q1 Rf1 Rw1]
for i = FDSelect
    % Step 2.1): define  Rf11 = Rf1(:,indd), Rf12 = Rf1(:,indf)
    indd = find(~SFDI(i,:)); indf = find(SFDI(i,:));
    QR.InputGroup.disturbances = p+mu+indd; % indices of Rf11
    QR.InputGroup.faults = p+mu+indf;       % indices of Rf12
    
    % Step 2.2): apply  Procedure EFD to {0,Rf11,Rf12,[Q1 Rf1 Rw1]}
    % to determine a least order Q1i such that Q1i*Rf11 = 0
    % QRfauxi contains the updated QRf: [Q1i*Rf12 Q1i*Q1 Q1i*Rf1 Q1i*Rw1]
    optionsi.HDesign = HDesign{i};
    if ~isempty(rdim)
       optionsi.rdim = rdim(i);
    end   
    try 
       [~,QRauxi,infoi] = efdsyn(QR,optionsi); 
    catch ME
       if strcmp(ME.message(1:5),'Empty') 
          error('Empty nullspace basis: the EFDIP is not solvable')
       else
          error(ME.message)
       end
    end
    QRi = QRauxi(:,'aux');   % extract [Q1i*Q1 Q1i*Rf1 Q1i*Rw1]
    
    % Step 2.3): 
    % extract Q{i} = Q1i*Q1 
    set(QRi,'InputGroup',struct('outputs',1:p,'controls',p+(1:mu)));
    Q{i} = QRi(:,1:p+mu);
    if nargout > 1
       % extract R{i} = [ Q1i*Rf1 Q1i*Rw1 ]
       set(QRi,'InputGroup',struct('faults',p+mu+(1:mf),'noise',p+mu+mf+(1:mw)));
       R{i} = QRi(:,p+mu+(1:mf+mw));
    end
    if nargout > 2
       info.tcond(i) = infoi.tcond;
       info.degs{i} = infoi.degs;
       info.HDesign{i} = infoi.HDesign;
    end       
end

% end EFDISYN
end

