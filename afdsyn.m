function [Q,R,info] = afdsyn(sysf,options)
%AFDSYN Approximate synthesis of fault detection filters
%  [Q,R,INFO] = AFDSYN(SYSF,OPTIONS) solves the approximate fault detection  
%  problem (AFDP) for a LTI system SYSF with additive faults.  
%  Two stable and proper filters Q and R are determined, where Q is the 
%  solution of the AFDP and R is the corresponding internal form. 
%
%  The continuous- or discrete-time system SYSF must be given in a standard
%  or descriptor state-space form, which corresponds to the input-output 
%  form  
%
%       y = Gu*u + Gd*d + Gf*f + Gw*w + Ga*aux,
%
%  with the Laplace- or Z-transformed plant outputs y, control inputs u, 
%  disturbance inputs d, fault inputs f, noise inputs w, and auxiliary 
%  inputs aux, and with Gu, Gd, Gf, Gw, and Ga the corresponding 
%  transfer-function matrices.
%  The inputs u, d, f, w and aux of SYSF correspond to five input groups 
%  named, respectively, {'controls','disturbances','faults','noise','aux'}.
%
%  The fault detection filter Q, determined in a standard state-space form, 
%  generates the residual signal r and corresponds to the 
%  input-output (implementation) form
%
%            r = Q*[ y ] = Qy*y + Qu*u .
%                  [ u ]
%
%  The inputs y and u of the resulting filter Q are grouped in two groups
%  {'outputs','controls'}, respectively, and the output group {'residuals'}
%  is defined for the residuals r. 
%
%  The filter R, determined in a standard state-space form, is the internal
%  form of Q, generates the residual signal r, and corresponds to the 
%  input-output form
%
%       r = Ru*u + Rd*d + Rf*f + Rw*w + Ra*aux ,
%
%  where 
%
%       [ Ru Rd Rf Rw Ra ] = Q*[ Gu Gd Gf Gw Ga ]. 
%                              [ I  0  0  0  0  ]
%
%  The solution of the AFDP ensures that Ru = 0, Rd = 0, Rf has all
%  its columns nonzero, while ||Rw||_inf < gamma 
%  (gamma is specified via the OPTIONS.gamma). 
%  The inputs f, w and aux of the resulting filter R are grouped 
%  in three input groups {'faults','noise','aux'}, respectively, and 
%  the output group {'residuals'} is defined for the residuals r. 

%  The resulting filters Q and R have, in general, the partitioned forms
%     Q = [ Q1 ] ,   R = [ R1 ] ,                      (1)
%         [ Q2 ]         [ R2 ]
%  where the filters Q1 and R1 with q1 outputs are a solution of the 
%  AFDP, while the filters Q2 and R2 with q2 outputs are a solution of an  
%  exact fault detection problem formulated for a reduced system obtained  
%  by decoupling the control and disturbance inputs from the residuals.  
%  The overall filters Q and R have observable state-space 
%  realizations (AQ,BQ,CQ,DQ) and (AQ,BR,CQ,DR), respectively, and thus  
%  share the observable pairs (AQ,CQ). 
%
%  The OPTIONS structure allows to specify various user options, as
%  follows:
%  OPTIONS.tol       - relative tolerance for rank computations
%                      (Default: internally determined value)
%  OPTIONS.tolmin    - absolute tolerance for observability tests
%                      (Default: internally determined value)
%  OPTIONS.FDTol     - threshold for fault detectability checks
%                      (Default: 0.0001))
%  OPTIONS.FDGainTol - threshold for strong fault detectability checks
%                      (Default: 0.01)
%  OPTIONS.rdim      - desired number q of residual outputs for Q and R
%                      (Default: [], in which case q = q1+q2, with q1 and 
%                      q2 chosen taking into account rw, the rank of the 
%                      transfer function matrix from the noise input to 
%                      the reduced system output, as follows:  
%                      if OPTIONS.HDesign is empty, then
%                         q1 = min(1,rw), if OPTIONS.minimal = true, or
%                         q1 = rw, if OPTIONS.minimal = false;
%                      if OPTIONS.HDesign is non-empty, then 
%                         q1 is the row dimension of the design 
%                         matrix H1 contained in OPTIONS.HDesign
%                      if OPTIONS.HDesign2 is empty, then
%                         q2 = 1-min(1,rw), if OPTIONS.minimal = true, or
%                         q2 = nvec-rw, where nvec is the dimension of the  
%                                left nullspace of the transfer function  
%                                matrix G1 := [Gu Gd;I 0], if 
%                                OPTIONS.minimal = false.)
%                      if OPTIONS.HDesign2 is non-empty, then 
%                         q2 is the row dimension of the design 
%                         matrix H2 contained in OPTIONS.HDesign2
%  OPTIONS.FDFreq    - vector of real frequency values for strong  
%                      detectability checks (Default: [])
%  OPTIONS.smarg     - stability margin for the poles of filters Q and R
%                      (Default: -sqrt(eps) for a continuous-time system 
%                               and 1-sqrt(eps) for a discrete-time system) 
%  OPTIONS.sdeg      - prescribed stability degree for the poles of the 
%                      filters Q and R
%                      (Default:  -0.05 for a continuous-time system and
%                                  0.95  for a discrete-time system) 
%  OPTIONS.poles     - complex vector containing a complex conjugate set  
%                      of desired poles (within the stability domain) 
%                      to be assigned for the filters Q and R (Default: []) 
%  OPTIONS.nullspace - specifies the proper nullspace basis option
%                      true  - use minimal proper basis (default); 
%                      false - use full-order observer based basis; 
%                              this option can be only used for a proper 
%                              system without disturbance inputs 
%  OPTIONS.simple    - option to employ a simple proper basis for synthesis 
%                      true  - use a simple basis; the orders of the  
%                              basis vectors are provided in INFO.deg
%                      false - no simple basis computed (default) 
%  OPTIONS.minimal   - option to perform least order filter synthesis 
%                      true  - perform least order synthesis (default)
%                      false - perform full order synthesis   
%  OPTIONS.exact     - option to perform exact filter synthesis 
%                      true  - perform exact synthesis 
%                      false - perform approximate synthesis (default)  
%  OPTIONS.tcond     - maximum alowed condition number of the employed 
%                      non-orthogonal transformations (Default: 1.e4).
%  OPTIONS.freq      - test frequency value to be employed to check the 
%                      full row rank admissibility condition
%                      (Default:[], i.e., a randomly generated frequency).
%  OPTIONS.HDesign   - design matrix H1, with full row rank q1,  to build  
%                      q1 linear combinations of the left nullspace basis 
%                      vectors of G1 := [ Gu Gd; I 0]. H1 is used in the 
%                      synthesis of the components Q1 and R1 in (1) 
%                      (Default: [])
%  OPTIONS.HDesign2  - design matrix H2, with full row rank q2, to build  
%                      q2 linear combinations of the left nullspace basis 
%                      vectors of G2 := [Gu Gd Gw; I 0 0]. H2 is used in  
%                      the synthesis of the components Q2 and R2 in (1). 
%                      (Default: []) 
%  OPTIONS.gamma     - upper bound on ||Rw||_inf  (Default: 1).
%  OPTIONS.epsreg    - regularization parameter   (Default: 0.1).
%  OPTIONS.sdegzer   - prescribed stability degree for zeros shifting  
%                      (Default:  -0.05 for a continuous-time system and
%                                  0.95  for a discrete-time system) 
%  OPTIONS.nonstd    - option to handle nonstandard optimization problems
%                      1 - use the quasi-co-outer-co-inner factorization
%                          (Default: [])
%                      2 - use the modified co-outer-co-inner factorization
%                          with the regularization parameter OPTIONS.epsreg
%                      3 - use the Wiener-Hopf type co-outer-co-inner 
%                          factorization
%                      4 - use the Wiener-Hopf type co-outer-co-inner 
%                          factorization with zero shifting of the 
%                          non-minimum phase factor using the
%                          stabilization parameter OPTIONS.sdegzer 
%                      5 - use the Wiener-Hopf type co-outer-co-inner 
%                          factorization with the regularization of the 
%                          non-minimum phase factor using the
%                          regularization parameter OPTIONS.epsreg
%
%  INFO is a structure containing additional information, as follows: 
%  INFO.tcond        - maximum of the condition numbers of the employed 
%                      non-orthogonal transformation matrices; a warning is 
%                      issued if INFO.tcond >= OPTIONS.tcond.
%  INFO.degs         - increasingly ordered degrees of a left minimal   
%                      polynomial nullspace basis of G1 := [ Gu Gd; I 0] 
%                      (also the left Kronecker indices of G1), if the 
%                      state-space realization of [Gu Gd] is minimal. 
%                      This information has been used in the case 
%                      OPTIONS.minimal = true to determine the least order
%                      components Q1 and R1 in (1).
%  INFO.degs2        - increasingly ordered degrees of a left minimal   
%                      polynomial nullspace basis of G2 := [Gu Gd Gw; I 0 0] 
%                      (also the left Kronecker indices of G2), if the 
%                      state-space realization of [Gu Gd Gw] is minimal.                      
%                      This information has been used in the case 
%                      OPTIONS.minimal = true to determine the least order
%                      components Q2 and R2 in (1).
%  INFO.S            - binary structure matrix corresponding to the  
%                      computed left nullspace basis of G1 = [ Gu Gd; I 0]
%                      used for the synthesis of the components Q1 and R1 
%                      in (1). 
%  INFO.S2           - binary structure matrix corresponding to the  
%                      computed left nullspace basis of 
%                      G2 = [ Gu Gd Gw; I 0 0] used for the synthesis of 
%                      the components Q2 and R2 in (1). 
%  INFO.HDesign      - the design matrix H1 employed for the synthesis of 
%                      the components Q1 and R1 in (1).
%  INFO.HDesign2     - the design matrix H2 employed for the synthesis of 
%                      the components Q2 and R2 in (1).
%  INFO.freq         - employed frequency FREQ to check the full row rank 
%                      admissibility condition.
%  INFO.gap          - achieved gap min_j ||Rf_j||_inf/||Rw||_inf.
%
%  See also EFDSYN, AFDISYN. 
%

%  Copyright 2018 A. Varga
%  Author:    A. Varga, 19-02-2018.
%  Revisions: A. Varga, 07-03-2018, 06-05-2018, 25-08-2018, 11-06-2019.
%
%  Method: The Procedure AFD from [1] is implemented, which is based 
%  on the synthesis method proposed in [2] and Remark 5.10 of [1]. 
%  The regularization approach based on the modified co-outer-co-inner 
%  factorization is discussed in [3] (see also Remark 5.8 of [1]). 
%
%  References:
%  [1] A. Varga
%      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%      Springer Verlag, 2017; sec. 5.3.
%  [2] A. Varga
%      General computational approach for optimal fault detection. 
%      Proc. IFAC Symposium SAFEPROCESS, Barcelona, Spain, pp. 107–112, 2009.
%  [3] K. Glover and A. Varga
%      On solving non-standard H-/H_2/inf fault detection problems. 
%      Proc. IEEE CDC, Orlando, FL, USA, pp. 891–896, 2011.

narginchk(1,2)
nargoutchk(0,3)

% check input system form
if ~isa(sysf,'ss')
   error('The input system SYSF must be an SS object')
end

if nargin < 2
    options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

discr = (sysf.Ts > 0);  % system type (continuous- or discrete-time)
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
      validateattributes(rdim, {'double'},{'integer','scalar','>=',0},'','OPTIONS.rdim')
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
strongFD = ~isempty(FDFreq);

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

% imposed design matrix H1 to form linear combinations of basis vectors
if isfield(options,'HDesign')
   HDesign1 = options.HDesign;
   if isempty(HDesign1)
      rdim1 = [];
   else 
      validateattributes(HDesign1, {'double'},{'2d'},'','OPTIONS.HDesign')
      rdim1 = size(HDesign1,1);
      if ~isempty(rdim) && rdim1 > rdim
         error('Row dimension of OPTIONS.HDesign must not exceed OPTIONS.rdim')
      end
      if rdim1 ~= rank(HDesign1)
         error('OPTIONS.HDesign must have full row rank')
      end
   end
else
   HDesign1 = []; rdim1 = [];
end
emptyHD1 = isempty(HDesign1);

% imposed design matrix H2 to form linear combinations of basis vectors
if isfield(options,'HDesign2')
   HDesign2 = options.HDesign2;
   if isempty(HDesign2)
      rdim2 = [];
   else
      validateattributes(HDesign2, {'double'},{'2d'},'','OPTIONS.HDesign2')
      rdim2 = size(HDesign2,1); 
      if ~isempty(rdim) && rdim2 > rdim
         error('Row dimension of OPTIONS.HDesign2 must not exceed OPTIONS.rdim')
      end
      if rdim2 ~= rank(HDesign2)
         error('OPTIONS.HDesign2 must have full row rank')
      end
   end
else
   HDesign2 = []; rdim2 = [];
end
emptyHD2 = isempty(HDesign2);

if ~isempty(rdim) && ~isempty(rdim1) && ~isempty(rdim2) && rdim1+rdim2 ~= rdim
    error('The sum of row dimensions of OPTIONS.HDesign and OPTIONS.HDesign2 must be equal to OPTIONS.rdim')
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
if isfield(sysf.InputGroup,'aux')
    % aux
    inpaux = sysf.InputGroup.aux;
    maux = length(inpaux);  
else
    inpaux = []; maux = 0;
end

m = mu+md+mf+mw+maux;       % total number of inputs
p = size(sysf,1);           % number of measurable outputs
lfreq = length(FDFreq);     % number of frequency values

if mf == 0 && minimal
   warning('Minimal synthesis option not feasible in the case of no faults')
   minimal = false;
end

% set default stability degree
if discr
   sdegdefault = 0.95;
else
   sdegdefault = -0.05;
end

% set options for nullspace computation
if strongFD
   % set options for nullspace computation with stabilization
   opts_glnull = struct('tol',tol,'m2',mf+mw+maux,'simple',simple,'sdeg',sdegdefault);
else 
   opts_glnull = struct('tol',tol,'m2',mf+mw+maux,'simple',simple);
end
% set options for LCF-based stabilization to be used for solvability checks
opts_glcf_default = struct('tol',tol,'tolmin',tolmin);
    
% Step 1): nullspace based reduction
%
if nullspace || md || rcond(sysf.e) < 1.e-7 
   % form [ Gu Gd Gf Gw Gaux; I 0 0 0 0] 
   syse = [sysf(:,[inpu inpd inpf inpw inpaux]); eye(mu,m)];
   %
   % compute a left nullspace basis Q = Q1 of G1 = [Gu Gd; I 0] = 0 and
   % obtain QR = [ Q R ], where R = [ Rf Rw Raux] = Q*[Gf Gw Ga;0 0 0]
   [QR,info1] = glnull(syse,opts_glnull); 
else
   % compute minimal basis as Q = Q1 = [ I -Gu] and set
   % QR = [ Q R ], where R = [ Gf Gw Ga ]
   QR = [ eye(p) dss(sysf.a,[-sysf.b(:,inpu) sysf.b(:,[inpf inpw inpaux])],...
          sysf.c,[-sysf.d(:,inpu) sysf.d(:,[inpf inpw inpaux])],sysf.e,sysf.Ts)];       
   if strongFD
      % perform stabilization if strong detectability has to be enforced
      QR = glcf(QR,opts_glcf_default);  
   end
   info1 = struct('degs',[],'tcond',1);
end


nvec = size(QR,1);          % number of basis vectors
% check solvability conditions
if nvec == 0,
   error('Empty nullspace basis: the AFDP is not solvable')
end
degs = info1.degs;          % degrees of a minimal polynomial basis
tcond1 = info1.tcond;       % condition number of employed transformations

indf = p+mu+(1:mf);           % input indices of Rf in QR
indw = p+mu+mf+(1:mw);        % input indices of Rw in QR
indaux = p+mu+mf+mw+(1:maux); % input indices of Ra in QR

% compute rank of Rw
rw = nrank(QR(:,indw),tol);


% set H1 for checking the solvability condition
if emptyHD1
   Htemp1 = eye(nvec);
else
   degs = [];
   [rdim1,nh] = size(HDesign1);
   if nh < nvec
      % pad with zeros: row rank is preserved
      Htemp1 = [ HDesign1 zeros(rdim1,nvec-nh) ];
   else
      % remove trailing columns if necessary: rank may drop
      Htemp1 = HDesign1(:,1:nvec);
      if nh > nvec && rdim1 ~= rank(Htemp1)
         error(['The leading ',num2str(rdim1),'x',num2str(nvec),...
                ' part of OPTIONS.HDesign must have full row rank'])
      end
      if rdim1 ~= nrank(Htemp1*QR(:,indw),tol)
         error('The addmissibility condition for OPTIONS.HDesign is not fulfilled')
      end
   end
end

% set H2 for checking the solvability condition
nvec2 = nvec-rw;
if nvec2
   if emptyHD2
      Htemp2 = eye(nvec2);
   else
      degs2 = [];
      [rdim2,nh] = size(HDesign2);
      if nh < nvec2
         % pad with zeros: row rank is preserved
         Htemp2 = [ HDesign2 zeros(rdim2,nvec2-nh) ];
      else
         % remove trailing columns if necessary: rank may drop
         Htemp2 = HDesign2(:,1:nvec2);
         if nh > nvec2 && rdim2 ~= rank(Htemp2)
            error(['The leading ',num2str(rdim2),'x',num2str(nvec2),...
                   ' part of OPTIONS.HDesign2 must have full row rank'])
         end
      end
   end
end

% adjust rdim, rdim1 and rdim2
if ~isempty(rdim), rdim = min(rdim,nvec); end
if ~isempty(rdim1), rdim1 = min(rdim1,rw); end
if ~isempty(rdim2), rdim2 = min(rdim2,nvec2); end

% setup the number of outputs of the filters Q1 and Q2 
if minimal
   %  least order design
   if isempty(rdim)
      % set default output dimensions 
      if emptyHD1 
         if rw 
            rdim1 = 1; 
         else
            rdim1 = 0; 
         end
      end
      if emptyHD2 
         if rw 
            rdim2 = 0; 
         else
            rdim2 = 1; 
         end
      end
   else
      % set output dimensions for a given rdim 
      if emptyHD1 && emptyHD2
         if rdim <= rw
            rdim1 = rdim; rdim2 = 0;
         else
            rdim1 = rw; rdim2 = rdim-rw;
         end
      elseif emptyHD1 
         rdim1 = rdim-rdim2;
      else
         rdim2 = rdim-rdim1;
      end
   end
else
   %  full order design
   if isempty(rdim) 
      if emptyHD1 
         rdim1 = rw;
      end
      if emptyHD2 
         rdim2 = nvec-rw;
      end
   else
      % set output dimensions for a given rdim 
      if emptyHD1 && emptyHD2
         if rdim <= rw
            rdim1 = rdim; rdim2 = 0;
         else
            rdim1 = rw; rdim2 = rdim-rw;
         end
      elseif emptyHD1 
         rdim1 = rdim-rdim2;
      else
         rdim2 = rdim-rdim1;
      end
   end
end

if rdim1
   % determine the structure matrix S1 underlying the synthesis of Q1
   if strongFD 
      S1 = fdisspec(Htemp1*QR(:,indf),FDGainTol,FDFreq);
   else
      S1 = fditspec(Htemp1*QR(:,indf),tol,FDTol);
   end
else
   if strongFD 
      S1 = false(0,mf,lfreq);
   else
      S1 = false(0,mf);
   end
end
    

if rdim2
   % the case rw < nvec
   % compute a left nullspace basis Qt such that Qt*Rw = 0 and
   % obtain QRt = [ Qt Rt ], where Rt = [ Q Rft Rwt Rat] = Qt*[Q1 Rf Rw Ra]
   opts_glnull.m2 = p+mu+mf+mw+maux;
   [QRt,infot] = glnull(QR(:,[indw 1:p+mu indf indw indaux]),opts_glnull);
   QR2 = QRt(:,nvec+[1:p+mu indf indw indaux]);
   if nvec2 ~= size(QR2,1);
      error('Something wrong: try to adapt the rank decision tolerance')
   end
   degs2 = infot.degs;
%
   % determine the structure matrix S2 underlying the synthesis of Q2
   if strongFD 
      S2 = fdisspec(Htemp2*QR2(:,indf),FDGainTol,FDFreq);
   else
      S2 = fditspec(Htemp2*QR2(:,indf),tol,FDTol);
   end
else
   if strongFD 
      S2 = false(0,mf,lfreq);
   else
      S2 = false(0,mf);
   end
end

% check solvability conditions
if nvec == 0,
   error('Empty nullspace basis: the AFDP is not solvable')
end


if isempty(indf)
   % handle the case of no faults as a normal case
   S1 = false(nvec,0);
   S2 = false(nvec2,0);
else
   S = [S1;S2];
   if strongFD 
      % check strong detectability conditions 
      for ii = 1:lfreq
          if ~all(max(S(:,:,ii),[],1))
             error('Strong detection of all faults not feasible')
          end  
      end
   else
      % check weak detectability conditions 
      if ~all(max(S,[],1))
         error('Detection of all faults not feasible')
      end   
   end
end

% determine offset value
if isempty(indf)
   foff = 0;
else
   foff = indf(1)-1;
end

if rdim1
   % synthesis of Q1 and R1
   QR1 = QR;
   indf1 = find(max(S1(:,:,1),[],1));
   options1 = options;
   QR1.InputGroup.faults = indf1+foff;
   QR1.InputGroup.noise = indw;
   options1.S = S1(:,indf1,:);
   options1.degs = degs; 
   options1.rdim = rdim1;
   [QR1,info1] = afdredsyn(QR1,options1);
   QR1.InputGroup.faults = indf;
   if nargout > 2
      info = info1;
      info.degs = degs;
      info.degs2 = [];
      info.S = S1;
      info.S2 = [];
      info.HDesign2 = [];
   end
else
   QR1 = gir(QR([],:));
   if nargout > 2
      info.tcond = 1;
      info.degs = [];
      info.S = S1;
      info.HDesign = [];
      info.gap = inf;
   end
end
if rdim2
   % synthesis of Q2 and R2
   indf2 = find(max(S2(:,:,1),[],1));
   QR2.InputGroup.faults = indf2+foff;
   QR2.InputGroup.noise = indw;
   QR2(:,indw) = 0*QR2(:,indw);
   options2 = options;
   options2.S = S2(:,indf2,:);
   options2.degs = degs2; 
   options2.rdim = rdim2;
   if emptyHD2
      options2.HDesign = [];
   else
      options2.HDesign = options.HDesign2;
   end
   [QR2,info2] = afdredsyn(QR2,options2);
   QR2.InputGroup.faults = indf;
   if nargout > 2
      info.tcond = max([tcond1;info.tcond;info2.tcond]);
      info.degs2 = degs2;
      info.S2 = S2;
      info.HDesign2 = info2.HDesign;
      info.gap = min(info.gap,info2.gap);
   end
else
   QR2 = gir(QR([],:));
   if nargout > 2
      info.tcond = 1;
      info.degs2 = [];
      info.S2 = S2;
      info.HDesign2 = [];
   end
end
QR = [QR1;QR2]; 

% set output variables
Q = QR(:,1:p+mu);
set(Q,'InputGroup',struct('outputs',1:p,'controls',p+(1:mu)));
set(Q,'OutputGroup',struct('residuals',1:size(Q,1)));

if nargout > 1
   R = QR(:,p+mu+1:end);
   set(R,'InputGroup',struct('faults',1:mf,'noise',mf+(1:mw),'aux',mf+mw+(1:maux)));
   set(R,'OutputGroup',struct('residuals',1:size(R,1)));
end

% end AFDSYN
end

function [sysredupd,info] = afdredsyn(sysred,options)
%AFDREDSYN Approximate synthesis of FD filters for a reduced system
%  [SYSREDUPD,INFO] = AFDREDSYN(SYSF,OPTIONS) determines Qred, the solution 
%  of the AFDP for a reduced system [Gf Gw] with Gf complete fault 
%  detectable and Gw a full-row rank transfer function matrix. 
%  The system SYSRED contains the state-space 
%  realization of [ Q Gf Gw Ga], where Q is the current FD filter, 
%  [Gf Gw] is the corresponding internal form (i.e., the reduced system), 
%  and Ga is an optional component. 
%  For the components Gf and Gw of SYSRED, the input groups 
%  {'faults'} and {'noise'} must be defined, respectively. 
%  The resulting system SYSREDUPD contains the state-space realization  
%  of the stable updated filters [ Qred*Q Qred*Gf Qred*Gw Qred*Ga] and
%  has the same input groups defined as for SYSRED. Note: Qred can be
%  explicitly determined, for example, by setting Ga = I. 
%
%  The solution Qred of the AFDP ensures that Qred*Gf has all
%  its columns nonzero, while ||Qred*Gw||_inf < gamma 
%  (gamma is specified via the OPTIONS.gamma). 
%
%  The OPTIONS structure allows to specify various user options, as
%  follows:
%  OPTIONS.tol       - relative tolerance for rank computations
%                      (Default: internally determined value)
%  OPTIONS.tolmin    - absolute tolerance for observability tests
%                      (Default: internally determined value)
%  OPTIONS.FDTol     - threshold for fault detectability checks
%                      (Default: 0.0001))
%  OPTIONS.FDGainTol - threshold for strong fault detectability checks
%                      (Default: 0.01)
%  OPTIONS.rdim      - desired number of residual outputs for Qred
%                      (Default: number of outputs of the reduced system)  
%  OPTIONS.FDFreq    - vector of real frequency values for strong  
%                      detectability checks (Default: [])
%  OPTIONS.smarg     - stability margin for the poles of filters Q and R
%                      (Default: -sqrt(eps) for a continuous-time system 
%                               and 1-sqrt(eps) for a discrete-time system) 
%  OPTIONS.sdeg      - prescribed stability degree for the poles of the 
%                      filters Q and R
%                      (Default:  -0.05 for a continuous-time system and
%                                  0.95  for a discrete-time system) 
%  OPTIONS.poles     - complex vector containing a complex conjugate set  
%                      of desired poles (within the stability domain) 
%                      to be assigned for the filters Q and R (Default: []) 
%  OPTIONS.simple    - option to employ a simple proper basis for synthesis 
%                      true  - use a simple basis; the orders of the  
%                              basis vectors are provided in INFO.deg
%                      false - no simple basis computed (default) 
%  OPTIONS.minimal   - option to perform a least order filter synthesis 
%                      true  - perform least order synthesis (default)
%                      false - perform full order synthesis   
%  OPTIONS.exact     - option to perform exact filter synthesis 
%                      true  - perform exact synthesis 
%                      false - perform approximate synthesis (default)  
%  OPTIONS.tcond     - maximum alowed condition number of the employed 
%                      non-orthogonal transformations (Default: 1.e4).
%  OPTIONS.HDesign   - full row rank design matrix H to build OPTIONS.rdim 
%                      linear combinations of the left nullspace basis 
%                      vectors (Default: [])
%  OPTIONS.gamma     - upper bound on ||Rw||_inf  (Default: 1).
%  OPTIONS.epsreg    - regularization parameter   (Default: 0.1).
%  OPTIONS.sdegzer   - prescribed stability degree for zeros shifting  
%                      (Default:  -0.05 for a continuous-time system and
%                                  0.95  for a discrete-time system) 
%  OPTIONS.nonstd    - option to handle nonstandard optimization problems
%                      1 - use the quasi-co-outer-co-inner factorization
%                          (Default: [])
%                      2 - use the modified co-outer-co-inner factorization
%                          with the regularization parameter OPTIONS.epsreg
%                      3 - use the Wiener-Hopf type co-outer-co-inner 
%                          factorization
%                      4 - use the Wiener-Hopf type co-outer-co-inner 
%                          factorization with zero shifting of the 
%                          non-minimum phase factor using the
%                          stabilization parameter OPTIONS.sdeg 
%                      5 - use the Wiener-Hopf type co-outer-co-inner 
%                          factorization with the regularization of the 
%                          non-minimum phase factor using the
%                          regularization parameter OPTIONS.epsreg
%  OPTIONS.degs      - increasingly ordered degrees of a left minimal   
%                      polynomial nullspace basis of G := [ Gu Gd; I 0] 
%                      (also the left Kronecker indices of G), if the 
%                      state-space realization of [Gu Gd ] is minimal                      
%  OPTIONS.S         - binary structure matrix corresponding to the  
%                      computed left nullspace basis.
%
%  INFO is a structure containing additional information, as follows: 
%  INFO.tcond        - maximum of the condition numbers of the employed 
%                      non-orthogonal transformation matrices; a warning is 
%                      issued if INFO.tcond >= OPTIONS.tcond.
%  INFO.HDesign      - the design matrix H employed for the synthesis of 
%                      the fault detection filter.
%  INFO.gap          - achieved gap min_j ||Rf_j||_inf/||Rw||_inf.

%  Copyright 2018 A. Varga
%  Author:    A. Varga, 04-05-2018.
%  Revisions: A. Varga, 25-08-2018,11-06-2019.
%
%  Method: The Procedure AFD from [1] is implemented, which is based 
%  on the synthesis method proposed in [2]. For the regularization approach 
%  based on the modified co-outer-co-inner factorization see [3]. For the
%  exact synthesis, the Procedure EFD from [1] is performed. 
%
%  References:
%  [1] A. Varga
%      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%      Springer Verlag, 2017; sec. 5.3.
%  [2] A. Varga
%      General computational approach for optimal fault detection. 
%      Proc. IFAC Symposium SAFEPROCESS, Barcelona, Spain, pp. 107–112, 2009.
%  [3] K. Glover and A. Varga
%      On solving non-standard H-/H_2/inf fault detection problems. 
%      Proc. IEEE CDC, Orlando, FL, USA, pp. 891–896, 2011.

narginchk(1,2)
nargoutchk(0,2)

% check input system form
if ~isa(sysred,'ss')
   error('The input system SYSF must be an SS object')
end
[nvec,m] = size(sysred);         
% nvec - number of reduced outputs
% m    - total number of inputs

% decode input information (only for the involved fault and noise inputs)
if isfield(sysred.InputGroup,'faults')
    % faults
    inpf = sysred.InputGroup.faults;
    mf = length(inpf);  
else
    inpf = []; mf = 0;
end
if mf == 0
%   error('No fault input group defined for SYSRED')
end

if isfield(sysred.InputGroup,'noise')
    % noise
    inpw = sysred.InputGroup.noise;
    mw = length(inpw);  
else
    inpw = []; mw = 0;
end


if nargin < 2
    options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

discr = (sysred.Ts > 0);  % system type (continuous- or discrete-time)
if discr
    smax = 1-sqrt(eps); smin = 0;
else
    smax = -sqrt(eps);  smin = -inf;
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
      validateattributes(rdim, {'double'},{'integer','scalar','>=',0},'','OPTIONS.rdim')
   end
else
   rdim = [];
end

% upper margin for condition  number of used transformations
if isfield(options,'tcond')
   tcond = options.tcond;
   validateattributes(tcond, {'double'},{'real','scalar','>=',1},'','OPTIONS.tcond') 
else
   tcond = 1.e4;
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
strongFD = ~isempty(FDFreq);

% desired stability degree
if isfield(options,'sdeg')
   sdeg = options.sdeg;
   if ~isempty(sdeg)
      validateattributes(sdeg, {'double'},{'real','scalar','<=',smax,'>=',smin},'','OPTIONS.sdeg') 
   end
else
   sdeg = [];
end
% set default stability degree
if discr
   sdegdefault = 0.95;
else
   sdegdefault = -0.05;
end

if isempty(sdeg)
   % set desired stability degree to default value
   sdeg = sdegdefault;
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


% option for simple basis
if isfield(options,'simple')
   simple = options.simple; 
   validateattributes(simple, {'logical'},{'binary'},'','OPTIONS.simple') 
else
   simple = false;
end

% information on nullspace basis vector degrees
if isfield(options,'degs')
   degs = options.degs; 
   if ~isempty(degs)
      validateattributes(degs, {'double'},{'vector','integer','>=',0},'','OPTIONS.degs') 
      if length(degs) ~= nvec
         error('Dimension of OPTIONS.degs must be equal to the number of reduced outputs')
      end
   end
else
   degs = [];
end

% structure matrix corresponding to employed basis
if isfield(options,'S')
   S = options.S; 
   if ~isempty(S)
      validateattributes(S, {'logical'},{'binary'},'','OPTIONS.S') 
      if size(S,2) ~= mf 
         error('Column dimension of OPTIONS.S must be equal to the number of faults')
      end
   end
else
   S = [];
end


% option for least order synthesis
if isfield(options,'minimal')
   minimal = options.minimal; 
   validateattributes(minimal, {'logical'},{'binary'},'','OPTIONS.minimal') 
else
    minimal = true;
end

% option for least order synthesis
if isfield(options,'exact')
   exact = options.exact; 
   validateattributes(exact, {'logical'},{'binary'},'','OPTIONS.exact') 
else
   exact = false;
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

% imposed design option to form linear combinations of basis vectors
if isfield(options,'HDesign')
   HDesign = options.HDesign;
   if ~isempty(HDesign)
      validateattributes(HDesign, {'double'},{'2d'},'','OPTIONS.HDesign')
      if ~isempty(rdim) && size(HDesign,1) ~= rdim
         error('Row dimension of OPTIONS.HDesign must be equal to OPTIONS.rdim')
      end
      if size(HDesign,1) ~= rank(HDesign)
         error('OPTIONS.HDesign must have full row rank')
      end
      if size(S,1) ~= size(HDesign,1) 
         error('Row dimensions of OPTIONS.S and OPTIONS.HDesign must be equal')
      end
   end
else
   HDesign = [];
end
emptyHD = isempty(HDesign);

% gamma, upper bound on ||Rw||_inf 
if isfield(options,'gamma')
   gamma = options.gamma;
   validateattributes(gamma, {'double'},{'real','scalar','>',0},'','OPTIONS.gamma') 
else
   gamma = 1;
end

% regularization parameter epsreg
if isfield(options,'epsreg')
   epsreg = options.epsreg;
   validateattributes(epsreg, {'double'},{'real','scalar','>',0},'','OPTIONS.epsreg') 
else
   epsreg = 0.1;
end

% desired stability degree for zero shifting
if isfield(options,'sdegzer')
   sdegzer = options.sdegzer;
   if ~isempty(sdegzer)
      validateattributes(sdegzer, {'double'},{'real','scalar','<=',smax,'>',smin},'','OPTIONS.sdegzer') 
   end
else
   sdegzer = [];
end
if isempty(sdegzer)
   % set desired stability degree of zeros to default value
   sdegzer = sdegdefault;
end

% option for handling nonstandard optimization problems
if isfield(options,'nonstd')
   jobio = options.nonstd;
   validateattributes(jobio, {'double'},{'integer','scalar','>=',1,'<=',5},'','OPTIONS.nonstd') 
else
   jobio = 1;
end


% set options for LCF-based stabilization to be used for admisibility checks
opts_glcf_default = struct('tol',tol,'tolmin',tolmin);
% set options for LCF-based stabilization to be used for final synthesis
opts_glcf = struct('tol',tol,'tolmin',tolmin, ...
                   'sdeg',sdeg,'smarg',smarg,'poles',poles);

lfreq = length(FDFreq);     % number of frequency values
tcond1 = 1;                 % condition number of employed transformations
nq = order(sysred);         % order of SYSRED

% set H for checking the solvability condition
if emptyHD
   Htemp = eye(nvec);
else
   degs = [];
   [rdim,nh] = size(HDesign);
   if nh < nvec
      % pad with zeros: row rank is preserved
      Htemp = [ HDesign zeros(rdim,nvec-nh) ];
   else
      % remove trailing columns if necessary: rank may drop
      Htemp = HDesign(:,1:nvec);
      if nh > nvec && rdim ~= rank(Htemp)
         error(['The leading ',num2str(rdim),'x',num2str(nvec),...
                ' part of OPTIONS.HDesign must have full row rank'])
      end
   end
end


% setup the number of filter outputs
if isempty(rdim)
   if emptyHD 
      rdim = 1;   
   else
      rdim = size(HDesign,1);
   end
else
   rdim = min(rdim,nvec); 
end


% no check of the solvability condition needed
if mw == 0
   rw = 0;
else
   rw = nrank(Htemp*sysred(:,inpw),tol);  
end

       
% Compute admissible Q2 to reduce the order of Q2*Q;  
% update Q <- Q2*Q, Rf = Q2*Gf 

sysredupd = sysred;
InpG = sysredupd.InputGroup; 

% reorder degs to correspond to the expected orders of basis vectors 
% corresponding to the actual order of outputs of QR 
if ~simple, degs = flip(degs); end
if rdim < nvec 
   % determine possible low order syntheses using i >= rmin basis vectors
   % and the corresponding expected orders    
   
   finish = false;    % set termination flag
   nout = rdim;       % initialize number of selected basis vectors
   if ~simple && minimal
      sysredupd = xperm(sysredupd,nq:-1:1);  % permute states to speedup glmcover1
   end
   if mw && ~exact
      rwgain = evalfr(Htemp*sysredupd(:,inpw),freq);
   else
      rwgain = Htemp(:,[]);  % set rwgain an empty matrix
   end
   itry = 1; 
   while ~finish     
       % choose nout basis vectors, which potentially lead to a least order
       % filter with rdim outputs:
       % basesel(i,:) contains the indices of candidate basis vectors;
       % ordsel(i)    contains the presumably achievable least orders
       [basesel,ordsel] = afdbasesel(S,rwgain,degs,rdim,nout,simple,tol);
       %
       % update the synthesis using the selections of candidate vector(s),
       % starting with the least (potentially) achievable order
       for i = 1:size(basesel,1);
           baseind = basesel(i,:); % indices of current basis selection
           if rdim == nout
               hbase = eye(rdim);
           else
               hbase = rand(rdim,nout); 
           end
           ip = [baseind, setdiff(1:nvec,baseind)]; 
           if simple
              if minimal
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
                 if rdim == nout
                    if emptyHD 
                       QRfwtest = modred(sysredupd(baseind,:),~noelim,'truncate');
                       h = Htemp(ip(1:rdim),:);
                    else
                       QRfwtest = gir(Htemp*sysredupd,tol);
                    end
                 else
                    % this case is possible only if HDesign is empty
                    % build rdim linear combinations of the first nout vectors 
                    QRfwtest = hbase*modred(sysredupd(baseind,:),~noelim,'truncate');
                    h = [ hbase zeros(rdim,nvec-nout) ]; 
                    h = h(:,ip);  % permute columns to match unpermuted QR 
                 end
              else
                 if rdim == nout
                    if emptyHD 
                       h = Htemp(ip(1:rdim),:);
                       QRfwtest = gir(sysredupd(baseind,:),tol,'finite');
                    else
                       QRfwtest = gir(Htemp*sysredupd,tol,'finite');
                    end
                 else
                    QRfwtest = gir(hbase*sysredupd(baseind,:),tol,'finite'); 
                    h = [ hbase zeros(rdim,nvec-nout) ]; 
                    h = h(:,ip);  % permute columns to match unpermuted QR 
                 end
              end
           else
              if minimal
                 % build output permutation vector for glmcover1  
                 if rdim == nout
                    if emptyHD 
                       [QRfwtest,info2] = glmcover1(sysredupd(ip,:),rdim,tol);
                       if ~isempty(ordsel) && (order(QRfwtest) ~= ordsel(i))
                          warning('AFDSYN: Expected reduced order not achieved')
                       end
                       h = Htemp(ip(1:rdim),:);
                    else
                       [QRfwtest,info2] = glmcover1([Htemp; eye(nvec)]*sysredupd(ip,:),rdim,tol);
                    end
                 else  
                    % this case is possible only if HDesign is empty
                    % build rdim linear combinations of the first nout vectors 
                    h = [ hbase zeros(rdim,nvec-nout) ]; 
                    [QRfwtest,info2] = glmcover1([h; eye(nvec)]*sysredupd(ip,:),rdim,tol);
                    h = h(:,ip);  % permute columns to match unpermuted QR 
                 end
              else
                 if rdim == nout
                    if emptyHD
                       h = Htemp(ip(1:rdim),:);
                       QRfwtest = gir(sysredupd(baseind,:),tol,'finite');
                    else
                       QRfwtest = gir(Htemp*sysredupd,tol,'finite');
                    end
                 else
                    QRfwtest = gir(hbase*sysredupd(baseind,:),tol,'finite'); 
                    h = [ hbase zeros(rdim,nvec-nout) ]; 
                    h = h(:,ip);  % permute columns to match unpermuted QR 
                 end
              end
           end
           % check admissibility of the current design; 
           if (rdim == nout && minimal) || rdim < nout
              % dismiss design if check fails
              if strongFD 
                 for ii = 1:lfreq
                     Stest = fdisspec(glcf(QRfwtest(:,inpf),opts_glcf_default),...
                                   FDGainTol,FDFreq(ii));
                     if ~all(max(Stest,[],1)), break, end
                 end
              else
                 Stest = fditspec(QRfwtest(:,inpf),tol,FDTol);
              end
              if rw 
                 if tol 
                     rwtest = rank(evalfr(QRfwtest(:,inpw),freq),tol);
                 else
                     rwtest = rank(evalfr(QRfwtest(:,inpw),freq));
                 end     
              else
                 rwtest = rdim;
              end
              % check complete fault detectability of the current design 
              % and full row rank condition  
              if all(max(Stest,[],1)) && rdim == rwtest
                 if ~simple && minimal
                    % adjust condition number of employed transformations
                    tcond1 = max([tcond1; info2.fnorm;info2.tcond]);
                    if tcond1 > tcond
                       disp(['AFDSYN: Possible loss of numerical stability',...
                            ' due to ill-conditioned transformations'])
                    end
                 end
                 sysredupd = QRfwtest;
                 finish = true;
                 break
              end
           else
              sysredupd = QRfwtest;
              finish = true;
              break
           end
       end
       nout = nout+1;
       if nout > nvec
          if itry > 5
             finish = true;
             warning('Fault detectability not achieved with the chosen number of residuals' )
          else
             itry = itry+1;
             nout = nout-1;
          end
       end
   end
   if emptyHD
      Htemp = h;
   end
else
   hbase = eye(rdim);
   if simple
      baseind = 1:rdim; 
   else
      baseind = 1;
   end
   h = eye(rdim);
   if ~emptyHD
      sysredupd = Htemp*sysredupd;
   else
      % use full minimum basis 
      Htemp = h;
   end
end

% compute M such that M*Q has a desired stability degree;  
% update Q <- M*Q and R <- M*R 
% this operation is performed only if rank is null or for exact synthesis
if rw == 0 || exact
   k = 1;
   if simple && isequal(hbase,eye(rdim)) && emptyHD 
      % exploit the block diagonal structure of basis matrices al and cl
      % to compute block-diagonal M
      [al,bl,cl,dl,el,Ts] = dssdata(sysredupd);
      for i = 1:length(baseind)
          blkord = degs(baseind(i));
          if blkord
             i1 = k:k+blkord-1; 
             QRfwi = glcf(dss(al(i1,i1),bl(i1,:),cl(i,i1),dl(i,:),el(i1,i1),Ts),opts_glcf);
             al(i1,i1) = QRfwi.a; bl(i1,:) = QRfwi.b;  cl(i,i1) = QRfwi.c;  
             dl(i,:) = QRfwi.d; 
             if isempty(QRfwi.e)
                el(i1,i1) = eye(blkord); 
             else
                el(i1,i1) = QRfwi.e;  
             end
             k = k+blkord;
          end
      end
      sysredupd = dss(al,bl,cl,dl,el,Ts);
   else
      sysredupd = glcf(sysredupd,opts_glcf);
   end
end
    
% finish if no noise input or if all noise inputs are decoupled or
% exact synthesis is performed
if mw == 0 || rw == 0 || exact
   % scale Rf to ensure unit minimum column gains
   if ~isempty(inpf)
      if strongFD && min(FDFreq) == 0
         % compute minimum DC gains  
         dcg = dcgain(sysredupd(:,inpf));  
         [y,indi] = max(abs(dcg),[],1);         % sort amplitudes of columns
         [scale,indj] = min(y);                 % select minimum amplitude
         sc = sign(dcg(indi(indj),indj))/scale; % adjust sign to be positive
      else
         % compute the minimum of H-inf norms of columns
         sc = 1/hinfminus(sysredupd(:,inpf));
      end
      sysredupd = sc*sysredupd;
   end

   % transform to standard state-space
   sysredupd = gss2ss(sysredupd,tol);
   set(sysredupd,'InputGroup',InpG)

   if nargout > 1
      % compute the achieved gap 
      if rw 
         beta = hinfminus(sysredupd(:,inpf),FDFreq);
         gap = beta/norm(sysredupd(:,inpw),inf);
      else
         gap = inf;
      end
      info = struct('tcond',tcond1,'HDesign',Htemp,'gap',gap); 
   end
   return
end

% determine the optimal factor Q3 to minimize the gap
% and update Q <- Q3*Q and R <- Q3*R 
  
% compute the extended quasi-co-outer-co-inner factorization  
% Rw = [Rwoe 0]*Rwi = Rwoe*Rwi1
[Rwi,Rwo,info1] = goifac(sysredupd(:,inpw),struct('tol',tol)); 
rw = size(Rwo,2);  % rank of Rw
nonstandard = info1.nfuz+info1.niuz > 0; 

% handle different cases
if nonstandard
   % non-standard case: zeros on the boundary of the stability domain 
   switch jobio
       case 1  
          % use quasi-co-outer-co-inner factorization
       case 2  
          % use modified co-outer-co-inner factorization
          [~,Rwo] = goifac([Rwo epsreg*eye(size(Rwo,1))],struct('tol',tol));
       case 3  
          % use Wiener-Hopf type co-outer-co-inner factorization
          % separate stable and strictly unstable zeros Rwo = Rwoe*Rwouz; 
          [Rwouz,Rwoe] = gcrange(Rwo,struct('tol',tol,'zeros','s-unstable'));
          Rwo = Rwoe*norm(Rwouz*Rwi(1:rw,:),inf);  
       case 4  
          % use modified Wiener-Hopf type co-outer-co-inner factorization
          % with zero shifting of the non-minimum phase factor
          % separate stable and strictly unstable zeros Rwo = Rwoe*Rwouz
          % and update Rwo <- Rwoe
          [Rwouz,Rwo] = gcrange(Rwo,struct('tol',tol,'zeros','s-unstable'));
          % set suitable bilinear transformation
          if sysred.Ts ~= 0
              sys1 = tf('z')/sdegzer;
          else
             s = tf('s'); sys1 = 1;
             if info1.nfuz
                sys1 = sys1*(s-sdegzer);
             end
             if info1.niuz
                sys1 = sys1/(1-sdegzer*s);
             end
          end
          % form shifted factor 
          Rwouz = gbilin(Rwouz,sys1); 
       case 5  
          % use modified Wiener-Hopf type co-outer-co-inner factorization
          % with regularization of the non-minimum phase factor
          % separate stable and strictly unstable zeros Rwo = Rwoe*Rwouz 
          % and update Rwo <- Rwoe
          [Rwouz,Rwo] = gcrange(Rwo,struct('tol',tol,'zeros','s-unstable'));
          rw = size(Rwouz,1);
          [~,Rwouz] = goifac([Rwouz epsreg*eye(rw)],struct('tol',tol));
   end
end

if rw == size(Rwo,1)
   % Q3 = inv(Rwo)
   % extract descriptor state-space data
   [aQR,bQR,cQR,dQR,eQR,Ts] = dssdata(sysredupd);
   % form [Rwo Q Rf Rw Raux] 
   RwoQR = dss(aQR,[Rwo.b bQR],cQR,[Rwo.d dQR],eQR,Ts);
   % form QR = inv(Rwo)*[ Q Rf Rw Raux] 
   sysredupd = grsol(RwoQR,m,struct('tol',tol));
   if nonstandard && (jobio == 4 || jobio == 5)
      sysredupd = gminreal(Rwouz\sysredupd,tol);
   end
else
   % regularization for non-invertible Rwo (this case should never occur)
   % Q3 = Rwoinv, where Rwoinv is a left inverse of Rwo
   % with stable spurious poles
   Rwoinv = glsol([Rwo;eye(rw)],rw,struct('tol',tol,'sdeg',sdeg)); 
   % form QR = Rwoinv*[ Q Rf Rw Raux] 
   sysredupd = gir(Rwoinv*sysredupd,tol,'finite');
   if nonstandard && (jobio == 4 || jobio == 5)
      sysredupd = gminreal(Rwouz\sysredupd,tol);
   end
end
if nonstandard && jobio == 1
   % perform stabilization 
   % determine Q4 such that Q <- Q4*Q and R <- Q4*R are stable
   sysredupd = glcf(sysredupd,opts_glcf);
end

% scale to enforce ||Rw||_inf = gamma
scale = gamma/norm(sysredupd(:,inpw),inf);
sysredupd = gss2ss(scale*sysredupd);
set(sysredupd,'InputGroup',InpG)

if nargout > 1
   % compute the achieved gap 
   beta = hinfminus(sysredupd(:,inpf),FDFreq);
   info = struct('tcond',tcond1,'HDesign',Htemp,'gap',beta/gamma); 
end

% end AFDREDSYN
end


