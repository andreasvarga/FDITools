function S = fdigenspec(sysf,options) 
%FDIGENSPEC  Generation of achievable FDI specifications.
%  S = FDIGENSPEC(SYSF,OPTION) determines, in the rows of the binary 
%  (logical) array S, all achievable FDI specifications for the
%  LTI state-space system SYSF with additive faults. 
%
%  The continuous- or discrete-time system SYSF must be given in a standard
%  or descriptor state-space form, which corresponds to the input-output 
%  form  
%
%       y = Gu*u + Gd*d + Gf*f ,
%
%  with the Laplace- or Z-transformed plant outputs y, control inputs u, 
%  disturbance inputs d, and fault inputs f, and with Gu, Gd, and Gf the 
%  corresponding transfer-function matrices.
%  The inputs u, d, and f of SYSF correspond to the three standard input 
%  groups named 'controls', 'disturbances', and 'faults', respectively. 
%  Each of these groups can be void. Any additionally defined input group
%  is ignored. 
%
%  If no standard input groups are explicitly defined, then an  
%  input-output representation, without control inputs, of the form
%
%                y = Gd*d + Gf*f 
%
%  is assumed, corresponding to a partition of SYSF as SYSF = [SYS1 SYS2], 
%  where Gd and Gf are the transfer-function matrices of SYS1 and SYS2, 
%  respectively. The dimension of the disturbance input d is specified
%  by the OPTIONS field OPTIONS.m1 (see below). If OPTIONS.m1 is specified,
%  then the above input-output form is assumed, even if the standard input
%  groups have been explicitly defined. 
%
%  The i-th row of S contains the i-th achievable specification,
%  obtainable by using a certain scalar output fault detection filter
%  with the input-output form
%
%                r = Q*[ y ] ,
%                      [ u ]
%
%  which ensures that in the corresponding internal form R = [Ru Rd Rf],
%  all control and disturbance inputs are decoupled from the residual r 
%  (i.e., Ru = 0 and Rd = 0) and the residual r is sensitive to a subset 
%  of fault inputs (i.e., Rf ~= 0).
%  The row S(i,:) is the structure matrix of Rf, whose element
%  S(i,j) = true if Rf(:,j) is nonzero and S(i,j) = false if Rf(:,j) = 0.
%  The check for nonzero Rf(:,j) is performed by using the function
%  FDITSPEC to evaluate the corresponding weak specifications.
%  If the frequency values for tests are provided (see below), then
%  |Rf(:,j)| must be larger than a threshold for all specified frequency
%  values. For this purpose, the function FDISSPEC is used to
%  evaluate the corresponding strong specifications.
%  The OPTIONS structure allows to specify various user options, 
%  as follows:
%  OPTIONS.tol       - tolerance for rank determinations
%                      (Default: internally computed)
%  OPTIONS.FDTol     - threshold for assessing weak specifications
%                      (see function FDITSPEC) (Default: 0.0001)
%  OPTIONS.FDGainTol - threshold for assessing strong specifications,
%                      i.e., threshold for nonzero frequency-responce
%                      gains for all frequency values specified in
%                      OPTIONS.FDFreq (see function FDISSPEC) 
%                      (Default: 0.01)
%  OPTIONS.m1        - M1, the number of inputs of SYS1 (Default: M1 = 0)
%                      If OPTIONS.m1 is specified, then SYSF is assumed to
%                      be partitioned as SYSF = [SYS1 SYS2] and the
%                      definitions of input groups are ignored. 
%  OPTIONS.FDFreq    - vector of real frequency values for strong  
%                      detectability checks (Default: [])
%  OPTIONS.sdeg      - prescribed stability degree for the
%                      poles of the generated filters (Default: [])
%
%  See also FDITSPEC and FDISSPEC.

%  Copyright 2016-2017 A. Varga
%  Author:      A. Varga, 25-11-2015.
%  Revision(s): A. Varga, 02-12-2016, 14-02-2017, 31-07-2017, 23-08-2017,
%                         24-10-2018.
%
%  Method: The Procedure GENSPEC from [1] is implemented. The nullspace 
%  method [2] is recursively employed to generate the complete set of 
%  achievable specifications. The method is also described in [3].
%
%  References:
%  [1] Varga A.
%      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%      Springer Verlag, 2017; sec. 5.4.
%  [2] Varga, A.:
%      On computing nullspace bases – a fault detection perspective. 
%      Proc. IFAC 2008 World Congress, Seoul, Korea, pages 6295–6300, 2008.
%  [3] Varga, A.:
%      On computing achievable fault signatures.
%      Proc. SAFEPROCESS'2009, Barcelona, Spain. 



% check input system form
if ~isa(sysf,'ss')
   error('The input system SYSF must be an SS object')
end

if nargin < 2
   options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

if isfield(options,'m1')
   % override input groups: apply inline procedure to G
   S = genspec_in(sysf,options);
   return
end

% decode input information to split G as G = [Gu Gd Gf]
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

m = mu+md+mf;

if m 
   % apply inline procedure to [ Gu Gd Gf; I 0 0 ] 
   options.m1 = mu+md;
   S = genspec_in([sysf(:,[inpu inpd inpf]); eye(mu,m)],options);
else
   % no input groups defined 
   % apply inline procedure to G
   options.m1 = 0;
   S = genspec_in(sysf,options);
end
   

% end FDIGENSPEC  
end
function S = genspec_in(sys,options) 
%GENSPEC_IN  Inline procedure to generate the achievable specifications.
%  S = GENSPEC_IN(SYS,OPTION) determines, in the rows of the binary
%  matrix S, all achievable fault detection specifications for the
%  partitioned state-space system SYS = [SYS1 SYS2], with the
%  corresponding input-output representation of the form,
%
%                y = Gd*d + Gf*f ,
%
%  where y is the output, d and f are the disturbance and fault
%  inputs, respectively, and Gd and Gf are the transfer-function
%  matrices of SYS1 and SYS2, respectively.
%  The i-th row of S contains the i-th achievable specification,
%  obtainable by using a certain scalar output fault detection filter
%  with the input-output form
%
%                r = Q*y ,
%
%  which ensures that all disturbance inputs d are decoupled
%  from the residual r (i.e., Q*Gd = 0) and the residual r is
%  sensitive to a subset of fault inputs (i.e., Q*Gf ~= 0).
%  The row S(i,:) is the structure matrix of Q*Gf, whose element
%  S(i,j) = true if Q*Gf(:,j) is nonzero and S(i,j) = false if Q*Gf(:,j) = 0.
%  The check for nonzero Q*Gf(:,j) is performed by using the function
%  FDITSPEC to evaluate the corresponding weak specifications.
%  If the frequency values for tests are provided (see below), then
%  |Q*Gf(:,j)| must be above a threshold for all specified frequency
%  values. For this purpose, the function FDISSPEC is used to
%  evaluate the corresponding strong specifications.
%  The OPTIONS structure allows to specify various user options, 
%  as follows:
%  OPTIONS.tol       - tolerance for rank determinations
%                      (Default: internally computed)
%  OPTIONS.FDTol     - threshold for assessing weak specifications
%                      (see function FDITSPEC) (Default: 0.0001)
%  OPTIONS.FDGainTol - threshold for assessing strong specifications,
%                      i.e., threshold for nonzero frequency-responce
%                      gains for all frequency values specified in
%                      OPTIONS.FDFreq (see function FDISSPEC) 
%                      (Default: 0.01)
%  OPTIONS.m1        - M1, the number of inputs of SYS1 (Default: M1 = 0)
%  OPTIONS.FDFreq    - vector of real frequency values for strong  
%                      detectability checks (Default: [])
%  OPTIONS.sdeg      - prescribed stability degree for the
%                      poles of the generated filters (Default: [])
%
%  See also FDITSPEC and FDISSPEC.

%  Copyright 2016 A. Varga
%  Author:      A. Varga, 25-11-2015.
%  Revision(s): A. Varga, 02-12-2016, 31-07-2017, 24-10-2018.
%
%  Method: The Procedure GENSPEC from [1] is implemented. The nullspace 
%  method [2] is recursively employed to generate the complete set of 
%  achievable specifications. The method is also described in [3].
%
%  References:
%  [1] Varga A.
%      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%      Springer Verlag, 2017; sec. 5.4.
%  [2] Varga, A.:
%      On computing nullspace bases – a fault detection perspective. 
%      Proc. IFAC 2008 World Congress, Seoul, Korea, pages 6295–6300, 2008.
%  [3] Varga, A.:
%      On computing achievable fault signatures.
%      Proc. SAFEPROCESS'2009, Barcelona, Spain. 


% set persistent variables used for recursive calls
persistent level tol sdeg smarg FDTol FDGainTol ksave Ts FDFreq nf

% determine input information (relevant for recursive usage)
if isfield(sys.InputGroup,'sys1')
   md = length(sys.InputGroup.sys1); 
else
   md = 0;
end
if isfield(sys.InputGroup,'sys2')
   mf = length(sys.InputGroup.sys2); 
else
   mf = 0;
end
       
% initialize recursion level
if md == 0 && mf == 0
   level = 0;   % iteration level is set to zero on first call
end

if level == 0
   % check input (only for level = 0) 
   if ~isa(sys,'ss')
      error('The input system SYS must be a state space LTI object')
   end   
   Ts = sys.ts;
   discr = (Ts > 0);
   m = size(sys,2);
   if discr
      smax = 1-sqrt(eps); smin = 0;
   else
      smax = -sqrt(eps);  smin = -inf;
   end
   
   % number of disturbance inputs
   if isfield(options,'m1')
      md = options.m1;
      validateattributes(md, {'double'},{'integer','scalar','>=',0,'<=',m},'','OPTIONS.m1') 
   else
      md = 0;
   end 
   mf = m-md;  % the number of fault inputs
 
   % finish if there are no faults
   if mf == 0
      S = [];
      return
   end

   % initialize OPTIONS if necessary
   if nargin == 1
      options = struct('tol',0);
   end
    
   % decode options (only for level = 0)
   
   % tolerance for rank determination
   if isfield(options,'tol')
      tol = max(0,options.tol);
      validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
   else
      tol = 0;
   end
   
   % tolerance for weak fault detectability checks
   if isfield(options,'FDTol')
      FDTol = options.FDTol;
      validateattributes(FDTol, {'double'},{'real','scalar','>=',0},'','OPTIONS.FDTol') 
   else
      FDTol = 0.0001;
   end
   
   % tolerance for strong fault detectability checks
   if isfield(options,'FDGainTol')
      FDGainTol = options.FDGainTol;
      validateattributes(FDGainTol, {'double'},{'real','scalar','>=',0},'','OPTIONS.FDGainTol') 
   else
      FDGainTol = 0.01;
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
   if ~isempty(sdeg) 
      if discr 
        smarg = 0.9;   % use a certain stability margin for poles
      else
        smarg = -0.05; % use a certain stability margin for poles
      end
   end
   if isempty(sdeg) 
      if discr 
        sdeg = 0.9;   % use a certain stability margin for poles
      else
        sdeg = -0.05; % use a certain stability margin for poles
      end
   end

   % frequencies for strong fault detectability checks
   if isfield(options,'FDFreq')
      FDFreq = options.FDFreq;
      if ~isempty(FDFreq)
         validateattributes(FDFreq, {'double'},{'real','vector','>=',0},'','OPTIONS.FDFreq')
      end
   else
      FDFreq = [];
   end
   nf = length(FDFreq);
   ksave = 1;  % parameter used to reduce the total number of iterations    
end

% compute SYSN*SYS2, where SYSN is a left nullspace of SYS1 
%    
if md > 0
   [redsys,~] = glnull(sys,struct('tol',tol,'m2',mf,'sdeg',sdeg));
   sys = redsys(:,end-mf+1:end);
end

% compute specification for the current level
if ~isempty(FDFreq) || ~isempty(sdeg) 
   % perform stabilization before evaluating current specification 
   % by computing a stable LCF
   warning('off','ACC:TEST')
   sysc = glcf(sys,struct('sdeg',sdeg,'smarg',smarg,'tolmin',tol));
   warning('on','ACC:TEST')
end

if isempty(FDFreq)
   % determine weak structure matrix
   if isempty(sdeg)  
      S = fditspec(sys,tol,FDTol); 
   else
      % use the stabilized system to evaluate S 
      S = fditspec(sysc,tol,FDTol); 
   end
else
   % determine strong structure matrix
   % enforce that strong specifications are a subset of weak specifications
   S = fditspec(sysc,tol,FDTol); S = [S; max(S,[],1)];
   S2 = fdisspec(sysc,FDGainTol,FDFreq); 
   for i = 1:nf
       %S = intersect(S,S2(:,:,i),'rows');
       S = intersect(S,[S2(:,:,i);max(S2(:,:,i),[],1)],'rows');
   end
end  
% add cummulated specification 
if size(S,1) > 1
   S = [S; max(S,[],1)];
end

% exit level if SYS has only one output, or if S = [] or S = 0
if size(sys,1) == 1  || isempty(S) || all(~max(S,[],1))   
   if level, level = level-1; end
   if size(S,1) > 1, S = unique(S,'rows'); S = S(any(S,2),:); end
   return 
end
%   loop over all fault inputs
for k1 = 1:mf
    if k1 < ksave, continue, end
    level = level+1;
    % redefine disturbance and fault inputs    
    sys.InputGroup.sys1 = 1; 
    sys.InputGroup.sys2 = 1:mf-1; 
    ksave = k1;
    s1 = genspec_in(sys(:,[k1, 1:k1-1,k1+1:mf])); 
    ksave = k1;
    %  add specifications from the current level
    S = [S; s1(:,1:k1-1),false(size(s1,1),1),s1(:,k1:mf-1)]; 
end
    
% eliminate duplicate or zero specifications
S = unique(S,'rows');
S = S(any(S,2),:); 
level = level-1;
    
% end GENSPEC_IN  
end



