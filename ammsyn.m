function [Q,R,info] = ammsyn(sysf,sysr,options)
%AMMSYN Approximate model matching based synthesis of FDI filters
%  [Q,R,INFO] = AMMSYN(SYSF,SYSR,OPTIONS) solves the approximate 
%  model-matching problem (AMMP) for a given LTI system SYSF with additive  
%  faults, to determine a stable fault detection and isolation filter Q,  
%  whose internal form approximates a given stable reference filter SYSR 
%  (possibly updated to enforce the stability of Q). 
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
%  The stable continuous- or discrete-time reference filter SYSR must be 
%  given in a standard or descriptor state-space form, which corresponds  
%  to the input-output form  
%
%       yr = Mru*u + Mrd*d + Mrf*f + Mrw*w ,
%
%  with the Laplace- or Z-transformed reference model outputs yr, control  
%  inputs u, disturbance inputs d, fault inputs f, and noise inputs w, and 
%  with Mru, Mrd, Mrf, and Mrw the corresponding transfer-function matrices.  
%  The inputs u, d, f, and w of SYSR correspond to four input groups 
%  named, respectively, {'controls','disturbances','faults','noise'}.
%
%  The stable filter Q, determined in a standard state-space form, 
%  corresponds to the input-output (implementation) form
%
%            r = Q*[ y ] = Qy*y + Q*u .
%                  [ u ]
%
%  The inputs y and u of the resulting filter Q are grouped in 
%  two input groups {'outputs','controls'}, respectively, and the output 
%  group {'residuals'} is defined for the residuals r. 
%
%  Let define
%              Ge = [ Gu Gd Gf Gw ],  Mr = [ Mru Mrd Mrf Mrw ] .                 
%                   [ I  0  0  0  ]
%
%  In the standard case, Ge has no zeros on the boundary of the 
%  stability domain, and the resulting stable filter is Q = Q0, where Q0 is
%  the optimal solution of the H-inf or H2-norm error minimization problem
%
%      gammaopt0 = ||Q0*Ge-M0*Mr|| = min,         (1)
%
%  where M0 is an updating factor. M0 = I in the case of emplyoing the 
%  H-inf norm, while in the case of employing the H2-norm, M0 = I for
%  a discrete-time system or, for a continuous-time system, 
%  M0 is determined a stable, diagonal, and invertible transfer function 
%  matrix, which ensures the existence of a finite H2-norm. 
%
%  In the non-standard case, Ge has zeros on the boundary of the 
%  stability domain, and the resulting optimal filter Q0, which solves 
%  the H-inf or H2-norm error minimization problem (1) is a possibly 
%  unstable or improper. A second updating factor M1 is determined, with   
%  the same properties as M0, which ensures that the computed stable and  
%  proper filter Q := M1*Q0 represents a suboptimal solution of an 
%  updated H-inf or H2-norm error minimization problem, for which the   
%  achieved suboptimal model-matching performance is
%
%      gammasub = ||Q*Ge-M1*M0*Mr|| .            (2)
%
%  As before, the optimal solution Qt of the updated H-inf or H2-norm error 
%  minimization problem 
%
%      gammaopt = ||Qt*Ge-M1*M0*Mr|| = min ,     (3)  
%
%  is possibly unstable or improper.
%
%  In both above cases, if G = [Gu Gd;I 0] has a non-empty left nullspace, 
%  then the nullspace method is employed to enforce the
%  exact decoupling of control and disturbance inputs from the residual
%  (i.e., Q*[Gu Gd;I 0] = 0 ). 
%
%  The stable filter R (also called the internal form of Q), represents the 
%  state-space realization of the filter R = Q*Ge := [ Ru Rd Rf Rw] and 
%  has the input groups {'controls','disturbances','faults','noise'}, where
%  only the input components corresponding to nonzero Ru, Rd, Rf or Rw 
%  are defined. 
%
%  If SYSR is a cell array which contains a bank of nb reference filters 
%  SYSR{1}, ..., SYSR{nb}, then each of the stable continuous- or 
%  discrete-time reference filter SYSR{i} must be given in a standard 
%  or descriptor state-space form, which corresponds to the input-output form  
%
%       yr_i = Mru_i*u + Mrd_i*d + Mrf_i*f + Mrw_i*w ,
%
%  with the Laplace- or Z-transformed reference model outputs yr_i, control  
%  inputs u, disturbance inputs d, fault inputs f, and noise inputs w, and 
%  with Mru_i, Mrd_i, Mrf_i, and Mrw_i the corresponding transfer-function 
%  matrices. The inputs u, d, f, and w of SYSR{i} correspond to four input 
%  groups named, respectively, {'controls','disturbances','faults','noise'}. 
%  Any of the input groups can be void, in which case, the corresponding 
%  transfer-function matrix is assumed to be zero. 
%
%  Q is a cell array which contains the stable filters Q{1}, ..., Q{nb},
%  where each Q{i} is determined in a standard state-space form and 
%  corresponds to the input-output (implementation) form
%
%            r_i = Q{i}*[ y ] = Qy_i*y + Qu_i*u .
%                       [ u ]
%
%  The inputs y and u of the resulting filter Q{i} are grouped in 
%  two input groups {'outputs','controls'}, respectively, and the output 
%  group {'residuals'} is defined for each residual component r_i. 
%
%  The filter Q{i} is the solution of an AMMP which involves
%
%      Ge = [ Gu Gd Gf Gw ],  Mr_i = [ Mru_i Mrd_i Mrf_i Mrw_i ] .                 
%           [ I  0  0  0  ]
%  
%  R is a cell array which contains the stable filters R{1}, ..., R{nb}, 
%  (also called the internal forms of Q{1}, ..., Q{nb}), and represents the 
%  state-space realization of the filter R{i} = Q{i}*Ge := [ Ru Rd Rf Rw].
%  Each R{i} has the input groups {'controls','disturbances','faults','noise'}, 
%  where only the input components corresponding to nonzero Ru, Rd, Rf or Rw 
%  are defined. 
%
%  The OPTIONS structure allows to specify various user options, as
%  follows:
%  OPTIONS.tol       - relative tolerance for rank computations
%                      (Default: internally determined value)
%  OPTIONS.tolmin    - absolute tolerance for observability tests
%                      (Default: internally determined value)
%  OPTIONS.reltol    - specifies the relative tolerance for the desired 
%                      accuracy of gamma-iteration. 
%                      (Default:  1.e-4)
%  OPTIONS.smarg     - prescribed stability margin for the poles of the 
%                      updating factor M
%                      (Default: -sqrt(eps) for a continuous-time system 
%                               and 1-sqrt(eps) for a discrete-time system) 
%  OPTIONS.sdeg      - prescribed stability degree for the poles of the 
%                      updating factor M
%                      (Default:  -0.05 for a continuous-time system and
%                                  0.95  for a discrete-time system) 
%  OPTIONS.poles     - specifies a complex conjugated set of desired poles 
%                      within the stability domain to be assigned for the 
%                      updating factor M (Default: []) 
%  OPTIONS.nullspace - specifies the proper nullspace basis option
%                      true  - use minimal proper basis (default); 
%                      false - use full-order observer based basis; 
%                              this option can be only used for a proper 
%                              system without disturbance inputs 
%  OPTIONS.simple    - option to employ a simple proper basis for synthesis 
%                      true  - use a simple basis; the orders of the  
%                              basis vectors are provided in INFO.deg
%                      false - no simple basis computed (default) 
%  OPTIONS.mindeg    - minimum degree solution 
%                      true  - determine, if possible, a  minimum order 
%                              stable solution  
%                      false - determine a particular stable solution which 
%                              has possibly non-minimal (default) 
%  OPTIONS.regmin    - specifies the regularization option with least order 
%                      left annihilator selection
%                      true  - perform least order selection (default)
%                      false - no least order selection to be performed  
%  OPTIONS.tcond     - maximum alowed condition number of the employed 
%                      non-orthogonal transformations (Default: 1.e4).
%  OPTIONS.normalize - specifies the normalization option for the diagonal 
%                      elements of the updating matrix M, as follows: 
%                      'gain'    - scale with the gains of the 
%                                  zero-pole-gain representation 
%                      'dcgain'  - scale with the DC-gains 
%                      'infnorm' - scale with the values of infinity-norm
%                                  (default)
%  OPTIONS.freq      - test frequency value to be employed to check the 
%                      left invertibility solvability condition
%                      (Default:[], i.e., a randomly generated frequency).
%  OPTIONS.HDesign   - full row rank design matrix H employed for the 
%                      synthesis of the filter Q (Default: []). 
%                      If SYSR contains nb filters, then separate design 
%                      matrices can be specified for each filter, as follows:
%                      OPTIONS.HDesign is an  nb-dimensional cell array, 
%                      where each OPTIONS.HDesign{i} is a full row rank 
%                      design matrix H_i employed for the synthesis of the 
%                      i-th fault detection filter Q{i} (Default: [])
%  OPTIONS.H2syn     - option to perform a H2-norm based synthesis 
%                      true  - perform a H2-norm based synthesis 
%                      false - perform a Hinf-norm based synthesis (default) 
%
%  INFO is a structure containing additional information, as follows: 
%  INFO.tcond        - maximum of the condition numbers of the employed 
%                      non-orthogonal transformation matrices; a warning is 
%                      issued if INFO.tcond >= OPTIONS.tcond.
%  INFO.degs         - increasingly ordered degrees of a left minimal   
%                      polynomial nullspace basis of G := [ Gu Gd; I 0] 
%                      (also the left Kronecker indices of G);
%  INFO.M            - employed updating matrix M in (1);
%                      the transfer-function matrix of M is stable, 
%                      diagonal and invertible.
%                      If SYSR contains nb filters, then INFO.M is an nb x 1  
%                      cell array, with INFO.M{i} containing the employed 
%                      stable, diagonal and invertible updating matrix used
%                      to determine a stable Q{i}. 
%  INFO.freq         - employed frequency FREQ to check left invertibility;
%  INFO.HDesign      - design matrix H employed for the synthesis of 
%                      the fault detection filter Q. 
%                      If SYSR contains nb filters, then INFO.HDesign is a 
%                      an nb x 1 cell array, with INFO.HDesign{i} 
%                      containing the design matrix H_i employed for the 
%                      synthesis of the filter Q{i}; 
%                      H = [] if no design matrix was involved. 
%  INFO.nonstandard  - logical value set as follows:
%                      true, for a non-standard problem (e.g., when 
%                            Ge(lambda) has zeros on the boundary of 
%                            the stability domain;
%                      false, for a standard problem (e.g., when 
%                            Ge(lambda) has no zeros on the boundary of 
%                            the stability domain;  
%                      If SYSR contains nb filters, then INFO.nonstandard 
%                      is a an nb x 1 logical vector, with 
%                      INFO.nonstandard(i) set as follows:
%                      true, for a non-standard i-th problem;
%                      false, for a standard i-th problem.
%  INFO.gammaopt0    - optimal performance value gammaopt0 for the original
%                      problem (1).
%                      If SYSR contains nb filters, then INFO.gammaopt0 
%                      is a an nb x 1 vector, with INFO.gammaopt0(i)
%                      containing the optimal performance value for the
%                      i-th original problem.  
%  INFO.gammaopt     - optimal performance value gammaopt for the updated
%                      problem (3);
%                      For the standard case, this value is gammaopt0, the 
%                      optimal approximation error in (1), while for a 
%                      non-standard problem this value is the optimal 
%                      approximation gammaopt error in (3) which 
%                      corresponds to an optimal filter Qt having possibly 
%                      poles on the boundary of stability domain. 
%                      If SYSR contains nb filters, then INFO.gammaopt 
%                      is a an nb x 1 vector, with INFO.gammaopt(i)
%                      containing the optimal performance value for the
%                      i-th updated problem if INFO.nonstandard(i) = true,
%                      or the i-th original problem if 
%                      INFO.nonstandard(i) = false. 
%  INFO.gammasub     - suboptimal performance value gammasub in (2). 
%                      For the standard case, this value is gammaopt0, the 
%                      optimal approximation error in (1).
%                      If SYSR contains nb filters, then INFO.gammasub 
%                      is a an nb x 1 vector, with INFO.gammasub(i)
%                      containing the suboptimal performance value for the
%                      i-th updated problem if INFO.nonstandard(i) = true,
%                      or is equalt to INFO.gammaopt0(i) if
%                      INFO.nonstandard(i) = false. 

%  Copyright 2018-2019 A. Varga
%  Author:    A. Varga, 21-02-2018.
%  Revisions: A. Varga, 06-03-2018, 06-07-2018, 13-08-2018, 24-08-2018,
%                       26-10-2018, 07-07-2019.
%
%  Method: The synthesis Procedure AMMS from [1] is implemented. The 
%  Procedure AMMS relies on the approximate model-matching synthesis method 
%  proposed in [2]. For more details on computational aspects see [3]. 
%
%  References:
%  [1] A. Varga, "Solving Fault Diagnosis Problems - Linear Synthesis 
%      Techniques", Springer Verlag, 2017; sec. 5.6.
%  [2] A. Varga. "Integrated computational algorithm for solving 
%      H_inf-optimal FDI problems". In Proc. of the IFAC World Congress, 
%      Milano, Italy, pp. 10187–10192, 2011.
%  [3] A. Varga. "Descriptor system techniques in solving H_2/H-Inf-optimal
%      fault detection and isolation problems". In L. T. Biegler,  
%      S. L. Campbell, and V. Mehrmann (Eds.), Control and Optimization 
%      with Differential-Algebraic Constraints, vol. 23 of Advances in 
%      Design and Control, pp. 105–125. SIAM, 2012.

narginchk(2,3)
nargoutchk(0,3)

% check input systems form
if ~isa(sysf,'ss')
   error('The input system SYSF must be an SS object')
end
Ts = sysf.Ts; 
discr = (Ts ~= 0);  % system type (continuous- or discrete-time)
if discr
    smax = 1-sqrt(eps); smin = 0;
    sdegnull = 0.95;
else
    smax = -sqrt(eps);  smin = -inf;
    sdegnull = -0.05;    
end

sysrss = isa(sysr,'ss');
if sysrss
   sysr = {sysr};
end

if isa(sysr,'cell')
   nb = length(sysr);
   rdim = ones(nb,1);
   validateattributes(sysr, {'cell'},{'vector'},'','SYSR')
   for i = 1:nb
       if ~isa(sysr{i},'ss')
          error('The input system SYSR must be cell array of LTI state space objects')
       end
       rdim(i) = size(sysr{i},1); % number of i-th residual outputs
   end
   for i = 1:nb
       if Ts ~= sysr{i}.Ts
          error('All components of the reference model SYSR must have the same sampling time')
       end
       if ~isempty(sysr{i}) && ~isstable(sysr{i})
          error('All components of the reference model SYSR must be stable systems')
       end
   end
else
   error('SYSR must be either a cell array or a LTI state space object')
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

if nargin < 3
    options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end
 
% decode options

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

% relative tolerance for the desired accuracy of gamma iteration
if isfield(options,'reltol')
   reltol = options.reltol;
   validateattributes(reltol, {'double'},{'real','scalar','>=',0},'','OPTIONS.reltol') 
else
   reltol = 1.e-4;
end

% upper margin for the condition number of used transformations
if isfield(options,'tcond')
   tcond = options.tcond;
   validateattributes(tcond, {'double'},{'real','scalar','>=',1},'','OPTIONS.tcond') 
else
   tcond = 1.e4;
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

% minimum degree solution option
if isfield(options,'mindeg')
   mindeg = options.mindeg;
   validateattributes(mindeg, {'logical'},{'binary'},'','OPTIONS.mindeg') 
else
   mindeg = false;
end

% option for least order basis selection
if isfield(options,'regmin')
   regmin = options.regmin; 
   validateattributes(regmin, {'logical'},{'binary'},'','OPTIONS.regmin') 
else
   regmin = true;
end

% option for normalization
if isfield(options,'normalize')
   normalize = options.normalize; 
   validateattributes(normalize, {'char'},{'nonempty'},'','OPTIONS.normalize') 
else
   normalize = 'infnorm';
end

switch normalize
    case 'gain'
       job = 0; 
    case 'dcgain'
       job = 1;
    case 'infnorm'
       job = 2;
    otherwise
       error('No such normalization option')
end

% frequency value for left invertibility checks
if isfield(options,'freq')
   freq = options.freq;
   if ~isempty(freq)
      validateattributes(freq, {'double'},{'scalar'},'','OPTIONS.freq')
   end
else
   freq = rand;
end

if isfield(options,'HDesign')
   HDesign = options.HDesign;
   if ~isempty(HDesign)
      if isa(HDesign,'double')
         if size(HDesign,1) ~= rank(HDesign)
            error('OPTIONS.HDesign must have full row rank')
         end
         [h{1:nb}] = deal(HDesign);
         HDesign = h;
      else        
         validateattributes(HDesign, {'cell'},{'vector','numel',nb},'','OPTIONS.HDesign')
      end
      for i = 1:nb
          rdimi = rdim(i);
          if ~isempty(HDesign{i})
              if size(HDesign{i},1) ~= rdimi
                 error(['Row dimension of OPTIONS.HDesign{',num2str(i),'} must be equal to ',num2str(rdimi)])
              end
              if rdimi <= mf 
                 if size(HDesign{i},1) ~= rank(HDesign{i})
                     error(['OPTIONS.HDesign{',num2str(i),'} must have full row rank'])
                 end
              else
                  warning(['Specified design matrix OPTIONS.HDesign{',num2str(i),'} ignored'])
                  HDesign{i} = [];
              end
          end
      end
   end
   emptyHD = all(cellfun('isempty',HDesign));
else
   HDesign = cell(1,nb);
   emptyHD = true;
end

% option to perform a H2-norm based synthesis
if isfield(options,'H2syn')
   H2syn = options.H2syn; 
   validateattributes(H2syn, {'logical'},{'binary'},'','OPTIONS.H2syn') 
   if H2syn
      normflag = 2;
   else
      normflag = inf;
   end
else
   H2syn = false;
   normflag = inf;
end

inpru = cell(nb,1);
inprd = cell(nb,1);
inprf = cell(nb,1);
inprw = cell(nb,1);

for i = 1:nb
   % decode input information for SYSR 
   if isfield(sysr{i}.InputGroup,'controls')
      % controls
      inpru{i} = sysr{i}.InputGroup.controls;
      if mu ~= length(inpru{i})  
         error('Incompatible control input groups between SYSF and SYSR')
      end
   end
   if isfield(sysr{i}.InputGroup,'disturbances')
      % disturbances
      inprd{i} = sysr{i}.InputGroup.disturbances;
      if md ~= length(inprd{i})  
         error('Incompatible disturbance input groups between SYSF and SYSR')
      end
   end
   if isfield(sysr{i}.InputGroup,'faults')
      % faults
      inprf{i} = sysr{i}.InputGroup.faults;
      if mf ~= length(inprf{i})  
         error('Incompatible fault input groups between SYSF and SYSR')
      end
   end
   if isfield(sysr{i}.InputGroup,'noise')
      % noise
      inprw{i} = sysr{i}.InputGroup.noise;
      if mw ~= length(inprw{i})  
         error('Incompatible noise input groups between SYSF and SYSR')
      end
   end
end

emptyru = all(cellfun('isempty',inpru));
emptyrd = all(cellfun('isempty',inprd));
emptyrw = all(cellfun('isempty',inprw));
strongfdi = false(nb,1); fd = false(nb,1); 

rd = gnrank(sysf(:,inpd),tol,freq);
SMr = fditspec(sysr,tol);

if emptyru &&  emptyrd &&  emptyrw 
   weakfdi = false(nb,1);  
   for ib = 1:nb
      rdimi = rdim(ib);
      rref = gnrank(sysr{ib}(:,inprf{ib}),tol,freq);
      if rdimi <= p-rd
         strongfdi(ib) = (rdimi == mf && rref == mf);
         if ~strongfdi(ib)
            fd(ib) = all(SMr(ib,:));
            weakfdi(ib) = ~fd(ib);
         end
      else
         weakfdi(ib) = true;
      end
   end
else
   weakfdi = true(nb,1); 
end

% set options for LCF-based stabilization to be used for final synthesis
opts_glcf = struct('tol',tol,'tolmin',tolmin, 'mininf',true,...
                   'mindeg',true,'sdeg',sdeg,'smarg',smarg,'poles',poles);
            
Q = cell(nb,1); R = cell(nb,1); M = cell(nb,1); Htemp = cell(nb,1);
nonstandard = false(nb,1); 
gopt0 = zeros(nb,1); gopt = zeros(nb,1); gsub = zeros(nb,1); 
ind = 1:nb;

   
if emptyru && emptyrd && emptyrw && any(~weakfdi)
   % address the case of void control, disturbance and noise inputs,
   % for square and invertible Mrf solve a strong FDI problem
        
   % Step 1) Nullspace based reduction
   %
   if nullspace || md || rcond(sysf.e) < 1.e-7 
      % form [ Gu Gd Gf Gw; I 0 0 0] 
      syse = [sysf(:,[inpu inpd inpf inpw]); eye(mu,m)];
      %
      % compute a left nullspace basis Q1 of G1 = [Gu Gd; I 0] = 0 and
      % obtain QR = [ Q1 R1 ], where R1 = [ Rf1 Rw1 ] = Q1*[Gf Gw;0 0]
      % set options for initial nullspace computation 
      opts_glnull = struct('tol',tol,'m2',mf+mw,'simple',simple,'sdeg',sdegnull);
      [QR,info1] = glnull(syse,opts_glnull);
      degs = info1.degs;      % degrees of a minimal polynomial basis
      tcond1 = info1.tcond;   % condition number of employed transformations
   else
      % compute minimal basis as Q1 = [ I -Gu] and set
      % QR = [ Q1 R1 ], where R1 = [ Rf1 Rw1 ] = [ Gf Gw ]
      QR = [ eye(p) dss(sysf.a,[-sysf.b(:,inpu) sysf.b(:,[inpf inpw])],...
           sysf.c,[-sysf.d(:,inpu) sysf.d(:,[inpf inpw])],sysf.e,sysf.Ts)];       
      degs = [];      
      tcond1 = 1;   
   end

   nvec = size(QR,1);       % number of basis vectors
   % check solvability conditions
   if nvec == 0,
      error('Empty nullspace basis: the nullspace-based approach to solve the AMMP is not applicable')
   end
   
   indf = p+mu+(1:mf);     % input indices of Rf1 in QR
   mfw = mf+mw;
   indR = p+mu+(1:mfw);     % input indices of R = [Rf1 Rw1] in QR
   
   

   Rftest = evalfr(QR(:,indf),freq);
   for ib = ind(~weakfdi)
     rdimi = rdim(ib);
     % set Htemp for checking the solvability condition
     if emptyHD
        Htemp{ib} = eye(nvec);
     else
        degs = [];
        nh = size(HDesign{ib},2);
        if nh < nvec
           % pad with zeros: row rank is preserved
           Htemp{ib} = [ HDesign{ib} zeros(rdimi,nvec-nh) ];
        else
           % remove trailing columns if necessary: rank may drop
           Htemp{ib} = HDesign{ib}(:,1:nvec);
           if nh > nvec && rdimi ~= rank(Htemp{ib})
                 error(['The leading ',num2str(mh),'x',num2str(nvec),...
                 ' part of OPTIONS.HDesign{',num2str(i),'} must have full row rank'])
           end
        end
     end   
%   
     if strongfdi(ib)
        % check strong isolability
        
        if tol > 0
           rf = rank(Htemp{ib}*Rftest,tol);
        else
           rf = rank(Htemp{ib}*Rftest);
        end
        if mf ~= rf
           warning('The system SYSF is not strongly isolable')
        end
     else
        % check complete detectability
        S = fditspec(Htemp{ib}*QR(:,indf),tol);
        if ~all(max(S,[],1))
           error('The system SYSF is not completely detectable')
        end   
     end
   
     if nvec > rdimi
        % 2) Determine Q2 such that Q2*Rf1 has maximal full row rank
        %    and Q2*Q1 has least McMillan degree
        if ~simple, degs = flip(degs); end
        nq = order(QR);
        finish = false;  % set termination flag
        nout = rdimi;  % set the number of basis vectors to be selected
        if ~simple && regmin && ~isempty(degs)
           % permute states to speedup glmcover1
           QR = xperm(QR,nq:-1:1);
        end
        if emptyHD
           h = eye(nvec);
        end
        while ~finish
           % choose nout basis vectors, which potentially lead to a
           % least order filter with rdim outputs:
           % basesel(i,:) contains the indices of candidate basis vectors;
           % ordsel(i)    contains the presumably achievable least orders
           [basesel,ordsel] = ammbasesel(Htemp{ib}*Rftest,degs,nout,simple,tol,strongfdi);
           % update the synthesis using the selections of candidate vector(s),
           % starting with the least (potentially) achievable order
           for i = 1:size(basesel,1);
               baseind = basesel(i,:); % indices of current basis selection
               if nout == mf
                  hbase = eye(mf);
               else
                  hbase = rand(mf,nout);
               end
               if simple
                  % handle simple basis
                  % here only the case rdim = nout can happen
                  ip = [baseind, setdiff(1:nvec,baseind)];
                  if regmin
                     if emptyHD
                        % select vectors and elliminate unobservable dynamics
                        noelim = false(nq,1);
                        ell = sum(degs(1:basesel(i,1)-1));
                        for jj = 1:nout
                            ellnext = sum(degs(1:baseind(jj)));
                            noelim(ell+1:ellnext) = true;
                            ell = ellnext;
                        end
                     end
                     if emptyHD
                        QRfwtest = modred(QR(baseind,:),~noelim,'truncate');
                        h = Htemp{ib}(ip(1:rdimi),:);
                     else
                        QRfwtest = gir(Htemp{ib}*QR(ip,:),tol);
                     end
                  else
                     if emptyHD
                        h = Htemp{ib}(ip(1:rdimi),:);
                        QRfwtest = gir(QR(baseind,:),tol,'finite');
                     else
                        QRfwtest = gir(Htemp{ib}*QR,tol,'finite');
                     end
                  end
               else
                  % handle minimal basis
                  % build output permutation vector for glmcover1
                  ip = [baseind, setdiff(1:nvec,baseind)];
                  if regmin
                     if rdimi == nout
                        if emptyHD
                           [QRfwtest,info2] = glmcover1(QR(ip,:),rdimi,tol);
                           if ~isempty(ordsel) && order(QRfwtest) ~= ordsel(i)
                              warning('AMMSYN: Expected reduced order not achieved')
                           end
                           h = Htemp{ib}(ip(1:rdimi),:);
                        else
                           [QRfwtest,info2] = glmcover1([Htemp{ib}; eye(nvec)]*QR(ip,:),rdimi,tol);
                        end
                     else
                        % the case rdimi < nout can only happen if no
                        % HDesign is explicitly provided
                        Htemp{ib} = blkdiag(hbase,eye(nvec-nout));
                        [QRfwtest,info2] = glmcover1([Htemp{ib}; eye(nvec)]*QR(ip,:),rdimi,tol);
                     end
                  else
                     % here only the case rdimi = nout can happen
                     if emptyHD
                        h = Htemp{ib}(ip(1:rdimi),:);
                        QRfwtest = gir(QR(baseind,:),tol,'finite');
                     else
                        QRfwtest = gir(Htemp{ib}*QR,tol,'finite');
                     end
                  end
               end
   
               % check admissibility of compressed Rf1;
               if ~simple && regmin
                  % dismiss design if the check fails
                  Rftest1 = evalfr(QRfwtest(:,indf),freq);
                  if strongfdi(ib) 
                     if (tol > 0 && rank(Rftest1,tol) == nout) || ...
                        (tol == 0 && rank(Rftest1) == nout)
                        finish = true;
                     end
                  else
                     beta = norm(Rftest1(:,1)); 
                     for j = 2:mf
                         beta  = min(beta,norm(Rftest1(:,j)));
                     end
                     if (tol > 0 && beta > tol) || ...
                        (tol == 0 && beta > nvec*mf*eps(norm(Rftest1)))
                        finish = true;
                     end
                  end
                  % here was possibly an error: fixed by moving outside
                  if finish    
                     % adjust condition number of employed transformations
                     tcond1 = max([tcond1; info2.fnorm;info2.tcond]);
                     if tcond1 > tcond
                        disp(['AMMSYN: Possible loss of numerical stability',...
                              ' due to ill-conditioned transformations'])
                     end
                     QR = QRfwtest;
                     break
                  end
               else
                  QR = QRfwtest;
                  finish = true;
                  break
               end
           end
           nout = nout+1;
           if nout > nvec && ~finish
              error('Something wrong: try perhaps with another test frequency')
           end
        end
        if emptyHD
           Htemp{ib} = h;
        end
     end
   
     % 3): determine Q3 = inv(Ro), from the extended quasi-co-outer-co-inner
     %     factorization R = [Ro 0]Ri;
     %     update Q <- Q3*Q and R <- Q3*R
     [Ri,Ro,info1] = goifac(QR(:,indR),struct('tol',tol));
     [lo,ro] = size(Ro);
   
     % detect nonstandard problem
     nonstandard(ib) = (info1.nfuz+info1.niuz > 0);
   
     if lo == ro && (order(Ro) == order(QR)) 
        % extract descriptor state-space data
        [ae,be,ce,de,ee,Ts] = dssdata(QR);
        % form [Ro Q R] (recall that Ro and QR share the same A, E and C matrices)
        % check concatanation formulas
        RoQR = dss(ae,[Ro.b be],ce,[Ro.d de],ee,Ts);
        % form QR = inv(Ro)*[Q R]
        QR = grsol(RoQR,p+mu+length(indR),struct('tol',tol));
     else
        % compute a left inverse of Ro with stable free poles
        Roinv = glsol([Ro;eye(ro)],ro,struct('tol',tol,'sdeg',sdeg)); 
        % update QR <- Roinv*QR 
        QR = gir(Roinv*QR,tol,'finite');
     end
   
     % form explicitly Rref = [ Mrf 0 ]
     Rref = [ sysr{ib}(:,inprf{ib}) zeros(rdimi,mw)];
   
     % compute F = [ F1 F2 ] =  Rref*Gi'
     F = Rref*Ri';
      
     % 4): Compute the solution Q4 of the LDP min||[F1-Q4 F2]||
     %     and update Q <- Q4*Q and R <- Q4*R
     if H2syn
        if tol
           told = tol;
        else
           told = 1.e4*eps;
        end
        if ~discr && norm(F.d(:,ro+1:end)) < told
           Mib = ss(eye(rdimi)); Mib.Ts = Ts; 
           F.d(:,ro+1:end) = zeros(size(F,1),mfw-ro);
        else
           % update Mr to achieve finite L2-norm
           Mib = ss(sdeg,1,-sdeg,0)*eye(rdimi);
           Rref = Mib*Rref;
           F = Rref*Ri';
        end
        % solve the H2-LDP min||[ F1-Q4 F2 ] ||_2
        [Qt,Qtu] = gsdec(F(:,1:ro),struct('tol',tol,'job','stable'));
        QR = Qt*QR;  %  update [ Q R ] <- Q4*[ Q R ]
        gopt0(ib) = norm(glcfid(gir([Qtu F(:,ro+1:end)],tol)),2);
     else
        % solve the H_inf-LDP min ||[ F1-Q4 F2 ] ||_inf
        opts = struct('tol',tol,'reltol',reltol);
        [Qt,gopt0(ib)] = glinfldp(F,mfw-ro,opts);
        QR = Qt*QR;  %  update [ Q R ] <- Q4*[ Q R ]
        Mib = ss(eye(rdimi)); Mib.Ts = Ts; 
     end
   
     if nonstandard(ib)
        % 5) compute diagonal Q5 such that Q5*Q has a desired stability degree;
        %    update Q <- Q5*Q and R <- Q5*R
        Qt = ss(zeros(size(QR)));
        for i = 1:rdimi
          % compute the LCF after removing unobservable eigenvalues
          [Qti,Mi] = glcf(gir(QR(i,:),tolmin,'obs'),opts_glcf);
          switch job
             % scale diagonal factors
             case 0  % scale with gain
                sc = get(zpk(Mi),'k');
             case 1  % scale with dcgain
                sc = dcgain(Mi);
             case 2  % scale with infinity norm
                sc = norm(Mi,inf)*sign(dcgain(Mi));
          end
          Qt(i,:) = Qti/sc; Mib(i,i) = Mib(i,i)*Mi/sc;
        end
        QR = Qt;
        % compute the optimal distance for the updated problem 
        if H2syn
           [~,Qtu] = gsdec(Mib*F(:,1:ro),struct('tol',tol,'job','stable'));
           gopt(ib) = norm(glcfid(gir([Qtu Mib*F(:,ro+1:end)],tol)),2);
        else
           [~,gopt(ib)] = glinfldp(Mib*F,mfw-ro,opts);
        end
     else
        gopt(ib) = gopt0(ib);
     end
     % transform to standard state-space
     QR = gss2ss(QR,tol);  
     Q{ib} = QR(:,1:p+mu); R{ib} = QR(:,p+mu+1:end);  
     M{ib} = Mib; 
   end
end
if any(weakfdi) 
   % apply the three-step procedure to solve the AMMP if SYSR has  
   % 'controls', 'disturbances' or 'noise' input groups
   degs = []; 
   finish = false;
   if emptyru && emptyrd  
      % perform nullspace-based synthesis if SYSR has no 'controls' and
      % 'disturbances' input groups
   
      % 1) Compute Q1, the left nullspace of [Gu Gd;I 0], and  
      %    R = [Rf1 Rw1] = Q1*[Gf Gw;0 0]
      %    QR contains Q1 in QR(:,1:p+mu) and R in QR(:,p+mu+1:end),
      %    where Rf1 is in QR(:,p+mu+1:p+mu+mf) and 
      %    Rw1 is in QR(:,p+mu+mf+1:end)
      mr = mf+mw;
      if nullspace || md || rcond(sysf.e) < 1.e-7 
         % form [ Gu Gd Gf Gw; I 0 0 0] 
         syse = [sysf(:,[inpu inpd inpf inpw]); eye(mu,m)];
         %
         % compute a left nullspace basis Q1 of G1 = [Gu Gd; I 0] = 0 and
         % obtain QR = [ Q1 R1 ], where R1 = [ Rf1 Rw1 ] = Q1*[Gf Gw;0 0]
         % set options for initial nullspace computation 
         opts_glnull = struct('tol',tol,'m2',mr,'simple',simple);
         [QR,info1] = glnull(syse,opts_glnull);
         tcond1 = info1.tcond;   % condition number of employed transformations
      else
         % compute minimal basis as Q1 = [ I -Gu] and set
         % QR = [ Q1 R1 ], where R1 = [ Rf1 Rw1 ] = [ Gf Gw ]
         QR = [ eye(p) dss(sysf.a,[-sysf.b(:,inpu) sysf.b(:,[inpf inpw])],...
               sysf.c,[-sysf.d(:,inpu) sysf.d(:,[inpf inpw])],sysf.e,sysf.Ts)];       
         tcond1 = 1;   
      end
      if size(QR,1)
          % non-empty nullspace
          finish = true;
          % separate implementation and internal forms
          Q1 = gir(QR(:,1:p+mu),tol); 
          R1 = QR(:,p+mu+(1:mr)); 
      end
   end
   if ~finish && emptyru 
      % 1) Compute Q1, the left nullspace of [Gu;I], and  
      %    R = [Rd1 Rf1 Rw1] = Q1*[Gd Gf Gw;0 0 0]
      %    QR contains Q1 in QR(:,1:p+mu) and R in QR(:,p+mu+1:end),
      %    where Rd1 is in QR(:,p+mu+1:p+mu+md), Rf1 is in 
      %    QR(:,p+mu+md+1:p+mu+md+mf) and Rw1 is in QR(:,p+mu+md+mf+1:end)
   
      % set options for nullspace computation
      mr = md+mf+mw;
      if nullspace || rcond(sysf.e) < 1.e-7 
         % form [ Gu Gd Gf Gw; I 0 0 0] 
         syse = [sysf(:,[inpu inpd inpf inpw]); eye(mu,m)];
         %
         % compute a left nullspace basis Q1 of G1 = [Gu; I] = 0 and
         % obtain QR = [ Q1 R1 ], where R1 = [ Rd1 Rf1 Rw1 ] = Q1*[Gd Gf Gw;0 0 0]
         % set options for initial nullspace computation 
         opts_glnull = struct('tol',tol,'m2',mr,'simple',simple);
         [QR,info1] = glnull(syse,opts_glnull);
         tcond1 = info1.tcond;   % condition number of employed transformations
      else
         % compute minimal basis as Q1 = [ I -Gu] and set
         % QR = [ Q1 R1 ], where R1 = [ Rf1 Rw1 ] = [ Gf Gw ]
         QR = [ eye(p) dss(sysf.a,[-sysf.b(:,inpu) sysf.b(:,[inpd inpf inpw])],...
               sysf.c,[-sysf.d(:,inpu) sysf.d(:,[inpd inpf inpw])],sysf.e,sysf.Ts)];       
         tcond1 = 1;   
      end
      if size(QR,1)
          % non-empty nullspace
          finish = true;
          % separate implementation and internal forms
          Q1 = gir(QR(:,1:p+mu),tol); 
          R1 = QR(:,p+mu+(1:mr)); 
      end
   end
   if ~finish && emptyrd 
      % 1) Compute Q1, the left nullspace of [Gd;0], and  
      %    R = [Ru1 Rf1 Rw1] = Q1*[Gu Gf Gw;I 0 0]
      %    QR contains Q1 in QR(:,1:p+mu) and R in QR(:,p+mu+1:end),
      %    where Ru1 is in QR(:,p+mu+1:p+mu+mu), Rf1 is in 
      %    QR(:,p+mu+mu+1:p+mu+mu+mf) and Rw1 is in QR(:,p+mu+mu+mf+1:end)
   
      % set options for nullspace computation
      mr = mu+mf+mw;
      % form [ Gd Gu Gf Gw; 0 I 0 0] 
      syse = [sysf(:,[inpd inpu inpf inpw]); zeros(mu,md) eye(mu,m-md)];
      %
      % compute a left nullspace basis Q1 of G1 = [Gd; 0] = 0 and obtain 
      % QR = [ Q1 R1 ], where R1 = [ Ru1 Rf1 Rw1 ] = Q1*[Ru Gf Gw;I 0 0]
      % set options for initial nullspace computation 
      opts_glnull = struct('tol',tol,'m2',mr,'simple',simple);
      [QR,info1] = glnull(syse,opts_glnull);
      tcond1 = info1.tcond;   % condition number of employed transformations
      if size(QR,1)
          % non-empty nullspace
          finish = true;
          % separate implementation and internal forms
          Q1 = gir(QR(:,1:p+mu),tol); 
          R1 = QR(:,p+mu+(1:mr)); 
      end
   end
   if ~finish
      % 1) Set Q1 = I and R = [Ru1 Rd1 Rf1 Rw1] = Q1*[Gu Gd Gf Gw;I 0 0 0]
      %    QR contains Q1 in QR(:,1:p+mu) and R in QR(:,p+mu+1:end),
      %    where Ru1 is in QR(:,p+mu+1:p+mu+mu), Rd1 is in 
      %    QR(:,p+mu+1:p+mu+mu+md), Rf1 is in QR(:,p+mu+mu+md+1:p+mu+mu+md+mf)
      %    and Rw1 is in QR(:,p+mu+mu+md+mf+1:end)
      mr = m;
      Q1 = ss(eye(p+mu)); Q1.Ts = Ts;
      R1 = [sysf(:,[inpu inpd inpf inpw]); eye(mu,m)];
      syse = R1;
      tcond1 = 1;
   end

   for ib = ind(weakfdi)
       rdimi = rdim(ib);
       % form explicitly Rref = [ Mru_i Mrd_i Mrf_i Mrw_i ] 
       rinp = zeros(0,m);
       if ~isempty(inpru{ib}) 
          rinp = [rinp;eye(mu,m)];
       end
       if ~isempty(inprd{ib}) 
          rinp = [rinp; zeros(md,mu) eye(md,m-mu)];
       end
       if ~isempty(inprf{ib}) 
          rinp = [rinp; zeros(mf,mu+md) eye(mf,m-mu-md)];
       end
       if ~isempty(inprw{ib}) 
          rinp = [rinp; zeros(mw,m-mw) eye(mw)];
       end 
       Rref = sysr{ib}(:,[inpru{ib} inprd{ib} inprf{ib} inprw{ib}])*rinp;
     
       % 2) Compute the approximate solution of Q2*Ge = Rref .
       opts_glasol = struct('tol',tol,'mindeg',mindeg,...
                        'reltol',reltol,'H2sol',H2syn,'sdeg',sdeg,'poles',poles);
       [Q2,info1] = glasol(R1,Rref(:,m-mr+1:end),opts_glasol); 
       Qib = gir(Q2*Q1,tol);
       gopt0(ib) = info1.mindist;
   
       % 3) compute diagonal M such that M*Q has a desired stability degree;  
       %    update Q <- M*Q
       Mib = ss(eye(rdimi)); Mib.Ts = Ts; 
       nonstandard(ib) = info1.nonstandard;
       if nonstandard(ib)
          %Qt = ss(zeros(size(Qib)));
          for i = 1:rdimi
             % compute the LCF after removing unobservable eigenvalues
             [Qti,Mi] = glcf(gir(Qib(i,:),tolmin,'obs'),opts_glcf);
             switch job
                % scale diagonal factors
                case 0  % scale with gain
                   sc = get(zpk(Mi),'k');
                case 1  % scale with dcgain
                   sc = dcgain(Mi);
                case 2  % scale with infinity norm
                sc = norm(Mi,inf)*sign(dcgain(Mi));
             end
             Qib(i,:) = Qti/sc; Mib(i,i) = Mib(i,i)*Mi/sc;
          end
          %Q = Qt;
          [~,info1] = glasol(R1,Mib*Rref(:,m-mr+1:end),opts_glasol); 
       end
       Q{ib} = gss2ss(Qib,tol); 
       % compute the internal form
       R{ib} = gir(Qib*syse(:,m-mr+1:end),tol); 
       gopt(ib) = info1.mindist; 
       M{ib} = Mib; 
   end
end

for ib = 1:nb
  rdimi = rdim(ib);
%   Qib = Q{ib};
%   % set output variables
%   set(Qib,'InputGroup',struct('outputs',1:p,'controls',p+(1:mu)));
%   set(Qib,'OutputGroup',struct('residuals',1:rdimi));
%   Q{ib} = Qib;
  Q{ib}.InputGroup.outputs = 1:p;
  Q{ib}.InputGroup.controls = p+(1:mu);
  Q{ib}.OutputGroup.residuals = 1:rdimi;

  if nargout > 1
     Rib = gss2ss(R{ib},tol); 
     set(Rib,'OutputGroup',struct('residuals',1:rdimi));
     ioff = 0;
     if ~isempty(inpru{ib})
        Rib.InputGroup.controls = 1:mu;
        ioff = ioff+mu;
     end
     if ~isempty(inprd{ib})
        Rib.InputGroup.disturbances = ioff+(1:md);
        ioff = ioff+md;
     end
     if ~isempty(inpf)
        Rib.InputGroup.faults = ioff+(1:mf);
        ioff = ioff+mf;
     end
     if ~isempty(inpw)
        Rib.InputGroup.noise = ioff+(1:mw);
     end
     R{ib} = Rib;
  end
end


if nargout > 2
   for ib = 1:nb
      if nonstandard(ib)
         gsub(ib) = fdimmperf(R{ib},M{ib}*sysr{ib},normflag); 
      else
         gsub(ib) = gopt(ib);
      end
   end
   if sysrss
       info = struct('tcond',tcond1,'degs',degs,'M',M{1},'freq',freq,...
                 'HDesign',Htemp{1},'gammaopt0',gopt0,'gammaopt',gopt, ...
                 'gammasub',gsub,'nonstandard',nonstandard); 
   else
       info = struct('tcond',tcond1,'degs',degs,'M',{M},'freq',freq,...
                 'HDesign',{Htemp},'gammaopt0',gopt0,'gammaopt',gopt, ...
                 'gammasub',gsub,'nonstandard',nonstandard); 
   end
end

if sysrss
   Q = Q{1}; R = R{1}; 
end


% end AMMSYN
end

