function [Q,R,info] = efdsyn(sysf,options)
%EFDSYN Exact synthesis of fault detection filters
%  [Q,R,INFO] = EFDSYN(SYSF,OPTIONS) solves the exact fault detection  
%  problem (EFDP) for a LTI system SYSF with additive faults.  
%  Two stable and proper filters Q and R are determined, where Q is the 
%  solution of the EFDP and R is the corresponding internal form. 
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
%  The solution of the EFDP ensures that Ru = 0, Rd = 0, and Rf has all
%  its columns nonzero. 
%  The inputs f, w and aux of the resulting filter R are grouped 
%  in three input groups {'faults','noise'}, respectively, and 
%  the output group {'residuals'} is defined for the residuals r. 
%
%  The resulting filters Q and R have observable state-space realizations
%  (AQ,BQ,CQ,DQ) and (AQ,BR,CQ,DR), respectively, and thus share the 
%  observable pairs (AQ,CQ). 
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
%                      (Default: [], in which case: 
%                                if OPTIONS.HDesign is empty, then
%                                   q = 1, if OPTIONS.minimal = true, or
%                                   q is the number of the nullspace basis  
%                                   vectors used for the initial synthesis, 
%                                   if OPTIONS.minimal = false;
%                                if OPTIONS.HDesign is non-empty, then 
%                                   q is the row dimension of the design 
%                                   matrix H contained in OPTIONS.HDesign
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
%  OPTIONS.minimal   - option to perform a least order filter synthesis 
%                      true  - perform least order synthesis (default)
%                      false - perform full order synthesis   
%  OPTIONS.tcond     - maximum alowed condition number of the employed 
%                      non-orthogonal transformations (Default: 1.e4).
%  OPTIONS.HDesign   - full row rank design matrix H to build OPTIONS.rdim 
%                      linear combinations of the left nullspace basis 
%                      vectors (Default: [])
%
%  INFO is a structure containing additional information, as follows: 
%  INFO.tcond        - maximum of the condition numbers of the employed 
%                      non-orthogonal transformation matrices; a warning is 
%                      issued if INFO.tcond >= OPTIONS.tcond.
%  INFO.degs         - increasingly ordered degrees of a left minimal   
%                      polynomial nullspace basis of G := [ Gu Gd; I 0] 
%                      (also the left Kronecker indices of G), if the 
%                      state-space realization of [Gu Gd ] is minimal                      
%  INFO.S            - binary structure matrix corresponding to the  
%                      computed left nullspace basis.
%  INFO.HDesign      - the design matrix H employed for the synthesis of 
%                      the fault detection filter.
%
%  See also EFDISYN. 
%

%  Copyright 2015-2018 A. Varga
%  Author:    A. Varga, 14-12-2015.
%  Revisions: A. Varga, 15-02-2016, 10-12-2016, 29-08-2017, 19-02-2018,
%                       05-06-2018, 07-06-2019. 
%
%  Method: The Procedure EFD from [1] is implemented. For more details on 
%  the least order synthesis of fault detection filters see [2].
%
%  References:
%  [1] Varga A.
%      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%      Springer Verlag, 2017; sec. 5.2.
%  [2] Varga, A.:
%      On computing least order fault detectors using rational 
%      nullspace bases. 
%      IFAC SAFEPROCESS'03 Symposium, Washington DC, USA, 2003.

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
   end
else
   HDesign = [];
end
emptyHD = isempty(HDesign);


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

% set options for nullspace computation
if strongFD
   % set options for nullspace computation with stabilization
   opts_glnull = struct('tol',tol,'m2',mf+mw+maux,'simple',simple,'sdeg',sdegdefault);
else 
   opts_glnull = struct('tol',tol,'m2',mf+mw+maux,'simple',simple);
end
% set options for LCF-based stabilization to be used for solvability checks
opts_glcf_default = struct('tol',tol,'tolmin',tolmin);
% set options for LCF-based stabilization to be used for final synthesis
opts_glcf = struct('tol',tol,'tolmin',tolmin, ...
                   'sdeg',sdeg,'smarg',smarg,'poles',poles);
    
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
if nvec == 0
   error('Empty nullspace basis: the EFDP is not solvable')
end
nq = order(QR);             % order of the minimal basis
degs = info1.degs;          % degrees of a minimal polynomial basis
tcond1 = info1.tcond;       % condition number of employed transformations

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
         error('The leading part of OPTIONS.HDesign must have full row rank')
      end
   end
end

indf = p+mu+(1:mf);   % input indices of Rf in QR
if isempty(indf)
   % handle the case of no faults as a normal case
   S = false(nvec,0);
else
   if strongFD 
      S = fdisspec(Htemp*QR(:,indf),FDGainTol,FDFreq);
      % check strong detectability conditions 
      for ii = 1:lfreq
          if ~all(max(S(:,:,ii),[],1))
             error('Strong detection of all faults not feasible')
          end  
      end
   else
      % check weak detectability conditions 
      S = fditspec(Htemp*QR(:,indf),tol,FDTol);
      if nvec == 0 || ~all(max(S,[],1))
         error('Detection of all faults not feasible')
      end   
   end
end

% setup the number of filter outputs
if minimal
   %  least order design
   if isempty(rdim)
      if emptyHD 
         rdim = 1;   
      else
         rdim = size(HDesign,1);
      end
   else
      rdim = min(rdim,nvec); 
   end
else
   %  full order design
   if isempty(rdim) 
      if emptyHD 
         rdim = nvec;
      else
         rdim = min(size(HDesign,1),nvec);
      end
   else
      if isempty(indf)
         if rdim < nvec && emptyHD
            warning(['rdim reset to ', num2str(nvec)])
            rdim = nvec;
         end
      else
         rdim = min(rdim,nvec); 
      end
   end
end
       
% Step 2): compute admissible Q2 to reduce the order of Q2*Q;  
% update Q <- Q2*Q, R <- Q2*R 

% reorder degs to correspond to the expected orders of basis vectors 
% corresponding to the actual order of outputs of QR 
if ~simple, degs = flip(degs); end
if rdim < nvec && ~isempty(indf)
   % determine possible low order syntheses using i >= rmin basis vectors
   % and the corresponding expected orders    
   
   finish = false;    % set termination flag
   nout = rdim;       % initialize number of selected basis vectors
   if ~simple && minimal
      QR = xperm(QR,nq:-1:1);  % permute states to speedup glmcover1
   end
   itry = 1; 
   while ~finish     
       % choose nout basis vectors, which potentially lead to a least order
       % filter with rdim outputs:
       % basesel(i,:) contains the indices of candidate basis vectors;
       % ordsel(i)    contains the presumably achievable least orders
       [basesel,ordsel] = efdbasesel(S,degs,rdim,nout,simple); 
       %
       % update the synthesis using the selections of candidate vector(s),
       % starting with the least (potentially) achievable order
       for i = 1:size(basesel,1)
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
                       QRfwtest = modred(QR(baseind,:),~noelim,'truncate');
                       h = Htemp(ip(1:rdim),:);
                    else
                       QRfwtest = gir(Htemp*QR,tol);
                    end
                 else
                    % this case is possible only if HDesign is empty
                    % build rdim linear combinations of the first nout vectors 
                    QRfwtest = hbase*modred(QR(baseind,:),~noelim,'truncate');
                    h = [ hbase zeros(rdim,nvec-nout) ]; 
                    h = h(:,ip);  % permute columns to match unpermuted QR 
                 end
              else
                 if rdim == nout
                    if emptyHD 
                       h = Htemp(ip(1:rdim),:);
                       QRfwtest = gir(QR(baseind,:),tol,'finite');
                    else
                       QRfwtest = gir(Htemp*QR,tol,'finite');
                    end
                 else
                    % this case is possible only if HDesign is empty
                    % build rdim linear combinations of the first nout vectors 
                    QRfwtest = gir(hbase*QR(baseind,:),tol,'finite'); 
                    h = [ hbase zeros(rdim,nvec-nout) ]; 
                    h = h(:,ip);  % permute columns to match unpermuted QR 
                 end
              end
           else
              if minimal
                 if rdim == nout
                    if emptyHD 
                       [QRfwtest,info2] = glmcover1(QR(ip,:),rdim,tol);
                       if ~isempty(ordsel) && (order(QRfwtest) ~= ordsel(i))
                          warning('EFDSYN: Expected reduced order not achieved')
                       end
                       h = Htemp(ip(1:rdim),:);
                    else
                       [QRfwtest,info2] = glmcover1([Htemp; eye(nvec)]*QR(ip,:),rdim,tol);
                    end
                 else  
                    % this case is possible only if HDesign is empty
                    % build rdim linear combinations of the first nout vectors 
                    h = [ hbase zeros(rdim,nvec-nout) ]; 
                    [QRfwtest,info2] = glmcover1([h; eye(nvec)]*QR(ip,:),rdim,tol);
                    h = h(:,ip);  % permute columns to match unpermuted QR 
                 end
              else
                 if rdim == nout
                    if emptyHD
                       h = Htemp(ip(1:rdim),:);
                       QRfwtest = gir(QR(baseind,:),tol,'finite');
                    else
                       QRfwtest = gir(Htemp*QR,tol,'finite');
                    end
                 else
                    QRfwtest = gir(hbase*QR(baseind,:),tol,'finite'); 
                    h = [ hbase zeros(rdim,nvec-nout) ]; 
                    h = h(:,ip);  % permute columns to match unpermuted QR 
                 end
              end
           end
           % check complete fault detectability of the current design; 
           if (rdim == nout && minimal) || rdim < nout
              % dismiss design if check fails
              if strongFD 
                 for ii = 1:lfreq
                     Stest = fdisspec(glcf(QRfwtest(:,indf),opts_glcf_default),...
                                   FDGainTol,FDFreq(ii));
                     if ~all(max(Stest,[],1)), break, end
                 end
              else
                 Stest = fditspec(QRfwtest(:,indf),tol,FDTol);
              end
              if all(max(Stest,[],1))
                 if ~simple && minimal
                    % adjust condition number of employed transformations
                    tcond1 = max([tcond1; info2.fnorm;info2.tcond]);
                    if tcond1 > tcond
                       disp(['EFDSYN: Possible loss of numerical stability',...
                            ' due to ill-conditioned transformations'])
                    end
%                     if ~emptyHD
%                        info.HDesign = Htemp;
%                     end
                 end
                 QR = QRfwtest;
                 finish = true;
                 break
              end
           else
              QR = QRfwtest;
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
      QR = Htemp*QR;
   else
      % use full minimum basis 
      Htemp = h;
   end
end

% Step 3): compute Q3 such that Q3*Q has a desired stability degree;  
% update Q <- Q3*Q, R <- Q3*R 
k = 1;
if simple && isequal(hbase,eye(rdim)) && emptyHD
    % exploit the block diagonal structure of basis matrices al and cl
    % to compute block-diagonal Q3
    [al,bl,cl,dl,el,Ts] = dssdata(QR);
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
    QR = dss(al,bl,cl,dl,el,Ts);
else
    QR = glcf(QR,opts_glcf);
end
    
% scale Rf to ensure unit minimum column gains
if ~isempty(indf)
   if strongFD && min(FDFreq) == 0
      % compute minimum DC gains  
      dcg = dcgain(QR(:,indf));  
      [y,indi] = max(abs(dcg),[],1);         % sort amplitudes of columns
      [scale,indj] = min(y);                 % select minimum amplitude
      sc = sign(dcg(indi(indj),indj))/scale; % adjust sign to be positive
   else
      % compute the minimum of H-inf norms of columns
      sc = 1/hinfminus(QR(:,indf));
   end
   QR = sc*QR;
end

% transform to standard state-space
QR = gss2ss(QR,tol);

% set output variables
Q = QR(:,1:p+mu);
set(Q,'InputGroup',struct('outputs',1:p,'controls',p+(1:mu)));
set(Q,'OutputGroup',struct('residuals',1:rdim));

if nargout > 1
   R = QR(:,p+mu+1:end);
   set(R,'InputGroup',struct('faults',1:mf,'noise',mf+(1:mw),'aux',mf+mw+(1:maux)));
   set(R,'OutputGroup',struct('residuals',1:rdim));
end

if nargout > 2
   info= struct('tcond',tcond1,'degs',info1.degs,'S',S,'HDesign',Htemp); 
end

% end EFDSYN
end

