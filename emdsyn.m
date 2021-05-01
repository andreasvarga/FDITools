function [Q,R,info]  = emdsyn(sysm,options) 
% EMDSYN Exact synthesis of model detection filters
%  [Q,R,INFO] = EMDSYN(SYSM,OPTIONS) solves the exact model detection  
%  problem (EMDP) for a multiple model SYSM with N component models.
%  Q is an N-vector of stable and proper filters, representing the 
%  solution of the EMDP. R is an NxN array of LTI systems, representing 
%  the corresponding internal forms. 
%
%  SYSM contains the multiple model specified either as an N-vector of 
%  LTI systems or a cell array with N elements, whose components are 
%  LTI systems with the same number of outputs and control inputs, and 
%  the same sampling time. 
%
%  The j-th component system is either a continuous- or discrete-time 
%  system in a standard or descriptor state-space form, which corresponds 
%  to the input-output form  
%
%      y_j = Gu_j*u + Gd_j*d_j + Gw_j*w_j ,
%
%  with the Laplace- or Z-transformed plant outputs y_j, control inputs u, 
%  disturbance inputs d_j, and noise inputs w_j, and with Gu_j, Gd_j 
%  and Gw_j the corresponding transfer-function matrices.
%  The inputs u, d_j, and w_j of each component system correspond to the  
%  input groups named, respectively, {'controls','disturbances','noise'}.
%
%  Q is an Nx1 cell array of filters, where the i-th filter Q{i},  
%  determined in a standard state-space form, generates the i-th residual  
%  signal r_i and corresponds to the input-output form
%
%            r_i = Q{i}*[ y ] = Qy_i*y + Qu_i*u ,
%                       [ u ]
%
%  where y and u are the measured ouput and control input of the 
%  underlying system.  The inputs y and u of each of the resulting filter 
%  Q{i} are grouped in two input groups {'outputs','controls'}, 
%  respectively, and the output group {'residuals'} is defined for 
%  the residuals r_i. 
%
%  R is an NxN cell array of filters, where the (i,j)-th filter R{i,j}, 
%  determined in a standard state-space form, is the internal form of Q{i}
%  acting on the j-th model, and corresponds to the input-output form
%
%       r_ij = Ru_ij*u + Rd_ij*d_j + Rw_ij*w_j ,
%
%  where 
%
%    R{i,j} := [Ru_ij Rd_ij Rw_ij] = Q{i}*[Gu_j Gd_j Gw_j]. 
%                                         [ I    0    0  ]
%
%  The solution of the EMDP ensures that Ru_ii = 0, Rd_ii = 0, and  
%  [Ru_ij Rd_ij] ~= 0 for i~= j.
%  The inputs u, d_j and w_j of each of the resulting filter R{i,j} are 
%  grouped in the input groups {'controls','disturbances','noise'}, 
%  respectively. 
%
%  The OPTIONS structure allows to specify various user options, as
%  follows:
%  OPTIONS.tol       - relative tolerance for rank computations
%                      (Default: internally determined value)
%  OPTIONS.tolmin    - absolute tolerance for observability tests
%                      (Default: internally determined value)
%  OPTIONS.MDTol     - threshold for model detectability checks
%                      (Default: 0.0001))
%  OPTIONS.MDGainTol - threshold for strong model detectability checks
%                      (Default: 0.01)
%  OPTIONS.rdim      - a vector q, whose i-th component q(i) specifies the 
%                      desired number of residual outputs for the filter
%                      Q{i}, or a scalar q, which specifies the   
%                      same desired number of residual outputs for all 
%                      filters Q{i} 
%                      (Default: [], in which case: 
%                                if OPTIONS.HDesign{i} is empty, then
%                                   q(i) = 1, if OPTIONS.minimal = true, or
%                                   q(i) is the dimension of the left  
%                                   nullspace which underlies the synthesis 
%                                   of Q{i}, if OPTIONS.minimal = false;
%                                if OPTIONS.HDesign{i} is non-empty, then 
%                                   q(i) is the row dimension of the i-th  
%                                   design matrix contained in 
%                                   OPTIONS.HDesign{i})
%  OPTIONS.MDFreq    - vector of nonnegative real frequency values for   
%                      strong model detectability checks (Default: [])
%  OPTIONS.emdtest   - specifies the option to perform extended model 
%                      detectability tests using both control and
%                      disturbance input channels 
%                      true  - use both control and disturbance input channels   
%                      false - use only the control channel (default)
%  OPTIONS.smarg     - prescribed stability margin for the poles of the 
%                      filters Q{i} 
%                      (Default: -sqrt(eps) for a continuous-time system 
%                               and 1-sqrt(eps) for a discrete-time system) 
%  OPTIONS.sdeg      - prescribed stability degree for the poles of the 
%                      filters Q{i} 
%                      (Default:  -0.05 for a continuous-time system and
%                                  0.95  for a discrete-time system) 
%  OPTIONS.poles     - specifies a complex conjugated set of desired poles 
%                      within the stability domain to be assigned for the 
%                      filters Q{i} (Default: []) 
%  OPTIONS.nullspace - specifies the nullspace option to employ in the 
%                      absence of disturbance inputs 
%                      true  - use minimal proper bases; 
%                      false - use simple observer based bases (default) 
%  OPTIONS.simple    - specifies the option to employ simple proper 
%                      left nullspace bases for synthesis 
%                      true  - use simple bases; the orders of the  
%                              basis vectors used to determine Q{i} 
%                              are provided in INFO.degs{i}
%                      false - no simple basis computed (default) 
%  OPTIONS.minimal   - specifies the option to perform least order 
%                      synthesis of the filters Q{i} and R{i,j}
%                      true  - perform least order synthesis (default)
%                      false - no least order synthesis performed  
%  OPTIONS.tcond     - maximum alowed condition number of the employed 
%                      non-orthogonal transformations (Default: 1.e4).
%  OPTIONS.MDSelect  - index set of desired filters to be designed
%                      (Default: 1:N)
%  OPTIONS.HDesign   - N-dimensional cell array; a nonempty full row 
%                      rank design matrix OPTIONS.HDesign{i} is a employed 
%                      for the synthesis of the i-th filter Q{i}; 
%                      if OPTIONS.rdim > 0, then each nonempty 
%                      OPTIONS.HDesign{i} must have at most OPTIONS.rdim 
%                      rows. 
%                      If OPTIONS.HDesign is specified as a full row rank 
%                      design matrix H, then an N-dimensional cell array
%                      with each OPTIONS.HDesign{i} = H is assumed.  
%                      (Default: []). 
%  OPTIONS.normalize - specifies the option to normalize the filters   
%                      Q{i} and R{i,j} such that the minimum gains of the
%                      off-diagonal elements of R{i,j} are equal to one;
%                      otherwise the standard normalization is  performed 
%                      to ensure equal gains for R{1,j} and R{j,1}
%                      true  - perform normalization to unit minimum gains
%                      false - perform standard normalization (default)  
%
%  INFO is a structure containing additional information, as follows: 
%  INFO.tcond        - maximum of the condition numbers of the employed 
%                      non-orthogonal transformation matrices; a warning is 
%                      issued if INFO.tcond >= OPTIONS.tcond.
%  INFO.degs         - N-dimensional cell array; INFO.degs{i} contains 
%                      the increasingly ordered left Kronecker indices of 
%                            G_i := [Gu_i Gd_i; I 0],  
%                      which are also the degrees of a left minimal   
%                      polynomial nullspace basis of G_i, and, if 
%                      OPTIONS.simple = true, also the orders of basis 
%                      vectors of the employed simple nullspace basis; 
%                      if OPTIONS.simple = false and Gd_i = 0, then 
%                      [I  -Gu_i] is used as a left nullspace basis and 
%                      INFO.degs{i} = []. 
%  INFO.MDperf       - NxN dimensional array containing the resulting  
%                      distance mapping performance, representing the peak  
%                      gains of the associated internal representations of 
%                      the model detection filters:
%                      if OPTIONS.MDFreq = [], then 
%                      INFO.MDperf(i,j) = ||[Ru_ij Rd_ij]||_inf if 
%                           OPTIONS.emdtest = true, or 
%                      INFO.MDperf(i,j) = ||Ru_ij||_inf, if 
%                            OPTIONS.emdtest = false;
%                      if OPTIONS.MDFreq ~= [], then INFO.MDperf(i,j) is 
%                      the maximum of the 2-norm gains of the frequency 
%                      responses of the selected input-output channels of 
%                      R{i,j} evaluated over all frequencies in 
%                      OPTIONS.MDFreq. Ideally, INFO.MDperf(i,j) should 
%                      map the distance between the i-th and j-th models 
%                      (e.g., as provided by the nugap metric).
%  INFO.HDesign      - N-dimensional cell array; INFO.HDesign{i} is the 
%                      design matrix H(i) employed for the synthesis of 
%                      the i-th filter.
%

%  Copyright 2017-2018 A. Varga
%  Author:      A. Varga, 28-08-2017.
%  Revision(s): A. Varga, 22-02-2018, 06-07-2018, 17-11-2018. 
%
%  Method:  The Procedure EMD from [1] is implemented, with the constraint
%  that all component model are stable and proper systems. For more details
%  on the least order synthesis of model detection filters see [2]. 
%
%  References:
%  [1] Varga A.
%      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%      Springer Verlag, 2017; sec. 6.2.
%  [2] Varga, A.:
%      Least order fault and model detection using multi-models. 
%      IEEE CDC'09, Shanghai, China, 2009.

narginchk(1,2)
nargoutchk(0,3)

if nargin < 2
    options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

if isa(sysm,'cell')
   N = length(sysm);
   validateattributes(sysm, {'cell'},{'vector'},'','SYSM')
   for i = 1:N
       if ~isa(sysm{i},'ss')
          error('The input system SYSM must be cell array of LTI state space objects')
       end
   end
   Ts = sysm{1}.Ts;
   for i = 2:N
       if Ts ~= sysm{i}.Ts
          error('All component models must have the same sampling time')
       end
   end
elseif isa(sysm,'ss')
   [~,~,N] = size(sysm);
   Ts = sysm.Ts; 
else
   error('SYSM must be either a cell array or a LTI state space object')
end

discr = (Ts ~= 0);
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
if isfield(options,'MDTol')
   MDTol = options.MDTol;
   validateattributes(MDTol, {'double'},{'real','scalar','>=',0},'','OPTIONS.MDTol') 
else
   MDTol = 0.0001;
end

% threshold for strong fault detectability checks
if isfield(options,'MDGainTol')
   MDGainTol =options.MDGainTol;
   validateattributes(MDGainTol, {'double'},{'real','scalar','>=',0},'','OPTIONS.MDGainTol') 
else
   MDGainTol = 0.01;
end

% desired number of filter outputs
if isfield(options,'rdim')
   rdim = options.rdim;
   if ~isempty(rdim)
      validateattributes(rdim, {'double'},{'integer','vector','>',0},'','OPTIONS.rdim')
      if length(rdim) == 1
         rdim = rdim*ones(1,N);
       else
         if N ~= length(rdim)
            error('The number of models must be equal to the dimension of the vector OPTIONS.rdim')
         end
       end
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
if isfield(options,'MDFreq')
   MDFreq = options.MDFreq;
   if ~isempty(MDFreq)
      validateattributes(MDFreq, {'double'},{'real','vector','>=',0},'','OPTIONS.MDFreq')
   end
else
   MDFreq = [];
end
strongMD = ~isempty(MDFreq);     % flag for strong model detection
lfreq = length(MDFreq);          % number of frequency values


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
   nullspace = false;
end


% option for simple basis
if isfield(options,'simple')
   simple = options.simple; 
   validateattributes(simple, {'logical'},{'binary'},'','OPTIONS.simple') 
else
   simple = false;
end
if simple
   % enforce nullspace option for simple basis
   nullspace = true; 
end 

% option for least order synthesis
if isfield(options,'minimal')
   minimal = options.minimal; 
   validateattributes(minimal, {'logical'},{'binary'},'','OPTIONS.minimal') 
else
   minimal = true;
end

% selected filters
if isfield(options,'MDSelect')
   MDSelect = options.MDSelect;
   validateattributes(MDSelect, {'double'},{'vector','integer','>=',1,'<=',N,'increasing'},'','OPTIONS.MDSelect')
else
   MDSelect = 1:N;
end

% imposed design matrix to form linear combinations of basis vectors
if isfield(options,'HDesign')
   HDesign = options.HDesign;
   if ~isempty(HDesign)
      if isa(HDesign,'double')
         if size(HDesign,1) ~= rank(HDesign)
            error('OPTIONS.HDesign must have full row rank')
         end
         [h{1:N}] = deal(HDesign);
         HDesign = h;
      else
         validateattributes(HDesign, {'cell'},{'vector','numel',N},'','OPTIONS.HDesign')
      end
      for i = 1:N
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
   HDesign = cell(1,N);
end

% option for extended model detectability test
if isfield(options,'emdtest')
   emdtest = options.emdtest; 
   validateattributes(emdtest, {'logical'},{'binary'},'','OPTIONS.emdtest') 
else
   emdtest = false;
end

% option for normalization
if isfield(options,'normalize')
   normalize = options.normalize; 
   validateattributes(normalize, {'logical'},{'binary'},'','OPTIONS.normalize') 
else
   normalize = false;
end
if normalize && length(MDSelect) < N 
   warning('Standard normalization disabled: using normalization to unit minimum gains')
   normalize = false; 
end
    
   



if isa(sysm,'ss') 
   [p,~,N]=size(sysm); 
   % decode input information
   if isfield(sysm.InputGroup,'controls')
      % controls
      tinpu = sysm.InputGroup.controls;
      mu = length(tinpu)*ones(N,1);  
   else
      tinpu = []; mu = zeros(N,1);
   end
   if isfield(sysm.InputGroup,'disturbances')
      % disturbances
      tinpd = sysm.InputGroup.disturbances;
      md = length(tinpd)*ones(N,1);  
   else
      tinpd = []; md = zeros(N,1);
   end
   if isfield(sysm.InputGroup,'noise')
      % noise
      tinpw = sysm.InputGroup.noise;
      mw = length(tinpw)*ones(N,1);  
   else
      tinpw = []; mw = zeros(N,1);
   end
   sys = cell(N,1);
   inpu = cell(N,1); inpd = cell(N,1); inpw = cell(N,1); 
   for i = 1:N
       sys{i} = sysm(:,:,i);
       inpu{i} = tinpu; 
       inpd{i} = tinpd; 
       inpw{i} = tinpw; 
   end
else
   N = length(sysm); 
   mu = zeros(N,1); md = zeros(N,1); mw = zeros(N,1);
   inpu = cell(N,1); inpd = cell(N,1); inpw = cell(N,1); 
   p = size(sysm{1},1); 
   Ts = sysm{1}.Ts;  sys = cell(N,1);
   for i = 1:N
       % decode input information
       if isfield(sysm{i}.InputGroup,'controls')
          % controls
          inpu{i} = sysm{i}.InputGroup.controls;
          mu(i) = length(inpu{i});  
       else
          inpu{i} = []; mu(i) = 0;
       end
       if isfield(sysm{i}.InputGroup,'disturbances')
          % disturbances
          inpd{i} = sysm{i}.InputGroup.disturbances;
          md{i} = length(inpd{i});  
       else
          inpd{i} = []; md(i) = 0;
       end
       if isfield(sysm{i}.InputGroup,'noise')
          % noise
          inpw{i} = sysm{i}.InputGroup.noise;
          mw{i} = length(inpw{i});  
       else
         inpw{i} = []; mw(i) = 0;
       end
       if i > 1 
          if p ~= size(sysm{i},1)
             error('All component models must have the same number of outputs')
          end
          if Ts ~= sysm{i}.Ts;
             error('All component models must have the same sampling time')
          end
       end
       sys{i} = sysm{i};
   end
   if max(mu) ~= min(mu)
      error('All component models must have the same number of control inputs')
   end       
end
 
% check stability of input model
for i = 1:N
    if ~isstable(sys{i})
       error('Input multiple model must be stable')
    end
end

% set options for nullspace computation
if minimal
   opt_glnull = struct('tol',tol,'simple',simple);
else
   opt_glnull = struct('tol',tol,'sdeg',sdeg,'poles',poles,'simple',simple);
end
% set options for LCF-based stabilization to be used for final synthesis
opts_glcf = struct('tol',tol,'tolmin',tolmin, ...
                   'sdeg',sdeg,'smarg',smarg,'poles',poles,'mininf',true);

Q = cell(N,1); 
if nargout > 1
   R = cell(N,N); 
end
if nargout > 2
   info.degs = cell(1,N); 
   info.HDesign = cell(1,N);
   info.tcond = ones(N,1); 
end
distinf = -ones(N,N); 
rdimtarget = rdim; 
for i = MDSelect
   % solve the exact MD problem for the i-th system (with noise inputs)
   % 
   % compute a minimal left nullspace basis Q1i for the i-th system
   if md(i) || nullspace
      % form Gi = [ Gui Gdi; I 0 ] for the i-th extended system                  
      syse = [ sys{i}(:,[inpu{i} inpd{i}]); eye(mu(i),mu(i)+md(i))]; 
      % compute a minimal left nullspace basis Q1i of Gi
      [qtemp,info1] = glnull(syse,opt_glnull);
   else
      % compute minimal basis as Q1i = [ I -Gui]
      qtemp = [eye(p) -sys{i}(:,'controls')];
      info1 = struct('degs',[],'tcond',1);
   end
   nvec = size(qtemp,1);       % number of basis vectors
   % check solvability conditions
   if nvec == 0,
      error('Empty nullspace basis: the EMDP is not solvable')
   end
   nq = order(qtemp);          % order of the minimal basis
   degs = info1.degs;          % degrees of a minimal polynomial basis
   if ~simple, degs = flip(degs); end
   tcond1 = info1.tcond;       % condition number of employed transformations
   
   % set H(i) for checking the solvability condition
   emptyHDi = isempty(HDesign{i});   
   if emptyHDi
      S = true(nvec,max(1,lfreq),N);
      Htemp = eye(nvec);
   else
      [mh,nh] = size(HDesign{i});
      S = true(mh,max(1,lfreq),N);
      degs = []; 
      if nh < nvec
         % pad with zeros: row rank is preserved
         Htemp = [ HDesign{i} zeros(mh,nvec-nh) ];
      else
         % remove trailing columns if necessary: rank may drop
         Htemp = HDesign{i}(:,1:nvec);
         if nh > nvec && mh ~= rank(Htemp)
            error(['The leading ',num2str(mh),'x',num2str(nvec),...
                   ' part of OPTIONS.HDesign{',num2str(i),'} must have full row rank'])
         end
      end
   end

   for j = 1:N
       % check j-th model detectability
       if i ~= j
          warning('off','all')
          temp = gir(Htemp*qtemp*[sys{j}(:,[inpu{j} inpd{j} inpw{j}]); eye(mu(j),mu(j)+md(j)+mw(j))],tol);
          warning('on','all')
          if strongMD
             if emdtest
                St = fdisspec(temp(:,1:mu(j)+md(j)),MDGainTol,MDFreq);
             else
                St = fdisspec(temp(:,1:mu(j)),MDGainTol,MDFreq);
             end
             Sj = true(size(St,1),lfreq); 
             for k = 1:lfreq
                 Sj(:,k) = any(St(:,:,k),2);
             end
             if ~any(max(Sj,[],1),2)
                if emdtest
                   error('Strong extended model detection not feasible')
                else
                   error('Strong model detection not feasible')
                end
             end 
          else
             if emdtest
                Sj = fditspec(temp(:,1:mu(j)+md(j)),tol,MDTol);
             else
                Sj = fditspec(temp(:,1:mu(j)),tol,MDTol);
             end
             if ~any(Sj(:))
                if emdtest
                   error('Extended model detection not feasible')
                else
                   error('Model detection not feasible')
                end
             end   
             Sj = max(Sj,[],2);
          end
          S(:,:,j) = Sj;
       end
   end
   
   % setup the number of filter outputs
   if minimal
      %  least order design
      if isempty(rdimtarget)
         if emptyHDi 
            rdim = 1;   
         else
            rdim = size(HDesign{i},1);
         end
      else
         rdim = min(rdimtarget(i),nvec); 
      end
   else
      %  full order design
      if isempty(rdimtarget)
         if emptyHDi 
            rdim = nvec;   
         else
            rdim = size(HDesign{i},1);
         end
      else
         rdim = min(rdimtarget(i),nvec); 
      end
   end

   if rdim < nvec
      % determine possible low order syntheses using i >= rmin basis vectors
      % and the corresponding expected orders 
   
      finish = false;    % set termination flag
      nout = rdim;       % initialize number of selected basis vectors
      if ~simple && minimal
         qtemp = xperm(qtemp,nq:-1:1);  % permute states to speedup glmcover1
      end
      
      % determine Q2i such that Q2i*Q1i is admisible and has least order      
      while ~finish
         % choose nout basis vectors, which potentially lead to a least order
         % filter with rdim outputs:
         % basesel(i,:) contains the indices of candidate basis vectors;
         % ordsel(i)    contains the presumably achievable least orders
         [basesel,ordsel] = emdbasesel(S,degs,rdim,nout,simple);
         %
         % update the synthesis using the selections of candidate vector(s),
         % starting with the least (potentially) achievable order
         for ii = 1:size(basesel,1);
             baseind = basesel(ii,:); % indices of current basis selection
             if rdim == nout
                hbase = eye(rdim);
             else
                hbase = rand(rdim,nout); 
             end
             ip = [baseind, setdiff(1:nvec,baseind)]; 
             if simple && ~isempty(degs)
                % select vectors and elliminate unobservable dynamics  
                if minimal
                   % determine Q2i such that Q2i*Q1i is admisible and 
                   % has least order
                   if emptyHDi 
                      noelim = false(nq,1); 
                      ell = sum(degs(1:basesel(ii,1)-1)); 
                      for jj = 1:nout 
                          ellnext = sum(degs(1:baseind(jj)));
                          noelim(ell+1:ellnext) = true;
                          ell = ellnext;
                      end
                   end
                   if rdim == nout
                      if emptyHDi 
                         qtest = modred(qtemp(baseind,:),~noelim,'truncate');
                         h = Htemp(ip(1:rdim),:);
                      else
                         qtest = gir(Htemp*qtemp,tol);
                      end
                   else
                      % this case is possible only if HDesign{i} is empty
                      % build rdim linear combinations of the first nout vectors 
                      qtest = hbase*modred(qtemp(baseind,:),~noelim,'truncate');
                      h = [ hbase zeros(rdim,nvec-nout) ]; 
                      h = h(:,ip);  % permute columns to match unpermuted qtemp 
                   end
                else 
                   % determine Q2i = Hi such that Q2i*Q1i is admisible 
                   if rdim == nout
                      if emptyHDi 
                         h = Htemp(ip(1:rdim),:);
                         qtest = gir(qtemp(baseind,:),tol,'finite');
                      else
                         qtest = gir(Htemp*qtemp,tol,'finite');
                      end
                   else
                      % this case is possible only if HDesign{i} is empty
                      % build rdim linear combinations of the first nout vectors 
                      qtest = gir(hbase*qtemp(baseind,:),tol,'finite'); 
                      h = [ hbase zeros(rdim,nvec-nout) ]; 
                      h = h(:,ip);  % permute columns to match unpermuted qtemp 
                   end
                end
             else
                if minimal
                   % determine Q2i such that Q2i*Q1i is admisible and 
                   % has least order
                   if rdim == nout
                      if emptyHDi 
                         [qtest,info2] = glmcover1(qtemp(ip,:),rdim,tol);
                         if ~isempty(ordsel) && (order(qtest) ~= ordsel(ii))
                            warning('EMDSYN: Expected reduced order not achieved')
                         end
                         h = Htemp(ip(1:rdim),:);
                      else
                         [qtest,info2] = glmcover1([Htemp;eye(nvec)]*qtemp(ip,:),rdim,tol);
                      end
                   else   
                      % this case is possible only if HDesign{i} is empty
                      % build rdim linear combinations of the first nout vectors 
                      h = [ hbase zeros(rdim,nvec-nout) ]; 
                      [qtest,info2] = glmcover1([h; eye(nvec)]*qtemp(ip,:),rdim,tol);
                      h = h(:,ip);  % permute columns to match unpermuted qtemp 
                   end
                else
                   % determine Q2i = Hi such that Q2i*Q1i is admisible 
                   if rdim == nout
                      if emptyHDi
                         h = Htemp(ip(1:rdim),:);
                         qtest = gir(qtemp(baseind,:),tol,'finite');
                      else
                         qtest = gir(Htemp*qtemp,tol,'finite');
                      end
                   else
                      qtest = gir(hbase*qtemp(baseind,:),tol,'finite'); 
                      h = [ hbase zeros(rdim,nvec-nout) ]; 
                      h = h(:,ip);  % permute columns to match unpermuted qtemp 
                   end
                end
             end
             
             % check model detectability of the current design; 
             if (rdim == nout && minimal) || rdim < nout
                notOK = false;
                for j = 1:N
                    % check i-th model detectability
                    if i ~= j
                       warning('off','all')
                       temp = gir(qtest*[sys{j}(:,[inpu{j} inpd{j} inpw{j}]);...
                                  eye(mu(j),mu(j)+md(j)+mw(j))],tol);
                       warning('on','all')
                       if strongMD
                          if emdtest
                             St = fdisspec(temp(:,1:mu(j)+md(j)),MDGainTol,MDFreq);
                          else
                             St = fdisspec(temp(:,1:mu(j)),MDGainTol,MDFreq);
                          end
                          Sj = true(size(St,1),lfreq); 
                          for k = 1:lfreq
                              Sj(:,k) = any(St(:,:,k),2);
                          end
                          notOK = ~any(max(Sj,[],1),2);
                       else
                          if emdtest
                             Sj = fditspec(temp(:,1:mu(j)+md(j)),tol,MDTol);
                          else
                             Sj = fditspec(temp(:,1:mu(j)),tol,MDTol);
                          end
                          notOK = ~any(Sj(:));
                       end 
                    end
                    if notOK, break, end
                end
                
                if ~notOK
                   if ~simple && minimal
                      % adjust condition number of employed transformations
                      tcond1 = max([tcond1; info2.fnorm;info2.tcond]);
                      if tcond1 > tcond
                         disp(['EMDSYN: Possible loss of numerical stability',...
                               ' due to ill-conditioned transformations'])
                      end
%                       if ~isempty(HDesign)
%                          info.HDesign{i} = Htemp;
%                       end
                   end
                   qtemp = qtest;
                   finish = true;
                   break
                end
             else
                qtemp = qtest;
                finish = true;
                break
             end
         end
         if ~finish
            if emptyHDi
               nout = nout+1;
               if nout > nvec
                  finish = true;
               end
            else
               qtemp = gir(Htemp*qtemp,tol);
               if minimal
                  warning('Model detection with least order not feasible for the given tolerances')
               end
               finish = true;
            end
         end
      end
      % set Q{i} = Q2i*Q1i 
      Q{i} = qtemp; 
      if emptyHDi
         Htemp = h;
      end
   else
      hbase = eye(rdim);
      baseind = 1:nvec;
      % set Q{i} = Q2i*Q1i 
      h = eye(rdim);
      if emptyHDi
         Q{i} = qtemp; 
         Htemp = h;
      else
         Q{i} = Htemp*qtemp; 
      end
   end
   
   % compute Q3i such that Q3i*Q{i} has a desired stability degree;  
   % update Q{i} <- Q3i*Qi{i} 
   k = 1;
   if simple && isequal(hbase,eye(rdim)) && emptyHDi
      % exploit the block diagonal structure of basis matrices al and cl
      % to compute block-diagonal Q3i
      [al,bl,cl,dl,el,Ts] = dssdata(Q{i});
      for ii = 1:length(baseind)
          blkord = degs(baseind(ii));
          if blkord
             i1 = k:k+blkord-1; 
             Qi = glcf(dss(al(i1,i1),bl(i1,:),cl(ii,i1),dl(ii,:),el(i1,i1),Ts),opts_glcf);
             al(i1,i1) = Qi.a; bl(i1,:) = Qi.b;  cl(ii,i1) = Qi.c;  
             dl(ii,:) = Qi.d; 
             if isempty(Qi.e)
                el(i1,i1) = eye(blkord); 
             else
                el(i1,i1) = Qi.e;  
             end
             k = k+blkord;
          end
      end
      Q{i} = dss(al,bl,cl,dl,el,Ts);
   else
      Q{i} = glcf(Q{i},opts_glcf);
   end

   for j = 1:N
       warning('off','all')
       temp = gir(Q{i}*[sys{j}(:,[inpu{j} inpd{j} inpw{j}]); eye(mu(j),mu(j)+md(j)+mw(j))],tol);
       warning('on','all')
       % compute gains
       if i ~= j
          if strongMD
             gains = 0;
             if emdtest
                for k = 1:lfreq
                    gains = max(gains,norm(evalfr(temp(:,1:mu(j)+md(j)),MDFreq(k))));
                end
             else
                for k = 1:lfreq
                    gains = max(gains,norm(evalfr(temp(:,1:mu(j)),MDFreq(k))));
                end
             end
             distinf(i,j) = gains; 
          else
             if emdtest
                distinf(i,j) = norm(temp(:,1:mu(j)+md(j)),inf);
             else
                distinf(i,j) = norm(temp(:,1:mu(j)),inf);
             end
          end
       else
          distinf(i,i) = 0;
       end
       if nargout > 1
          R{i,j} = temp;
       end
   end
   
   if normalize
      % scale to unit minimum gains
      scale = 1/min(distinf(i,[1:i-1 i+1:N]));
   else
      % scale to ensure distinf(i,1) = distinf(1,i)
      if i == 1
         scale = 1/min(distinf(1,2:N));
      else
         scale = distinf(1,i)/distinf(i,1);
      end
   end
   distinf(i,:) = scale*distinf(i,:);
      
   Q{i} = scale*gss2ss(Q{i});
   Q{i}.InputGroup = struct('outputs',1:p,'controls',p+(1:mu));
   Q{i}.OutputGroup = struct('residuals',1:rdim);
   if nargout > 1
      for j = 1:N
          R{i,j} = scale*gss2ss(R{i,j});
          R{i,j}.OutputGroup = struct('residuals',1:rdim);
      end
   end
   
   if nargout > 2
      info.degs{i} = info1.degs;  
      info.tcond(i) = tcond1;
      if emptyHDi
         info.HDesign{i} = h;
      else
         info.HDesign{i} = Htemp;
      end
   end
end

if nargout > 2
   info.MDperf = distinf;
end

% end EMDSYN
