function [Q,R,info] = emmsyn(sysf,sysr,options)
%EMMSYN Exact model matching based synthesis of FDI filters
%  [Q,R,INFO] = EMMSYN(SYSF,SYSR,OPTIONS) solves the exact model-matching 
%  problem (EMMP) for a given LTI system SYSF with additive faults, 
%  to determine a stable fault detection and isolation filter Q
%  such that its internal form R is equal to a given stable reference 
%  filter SYSR (possibly updated to enforce the stability of Q). 
%  If SYSR is specified as a bank of nb reference filters 
%  SYSR{1}, ..., SYSR{nb}, then the resulting Q contains a bank of filters 
%  Q{1}, ..., Q{nb}, whose internal forms  R{1}, ..., R{nb} are equal to
%  the corresponding specified reference filters.
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
%  If SYSR contains a stable continuous- or discrete-time reference filter,
%  then SYSR must be given in a standard or descriptor state-space form, 
%  which corresponds to the input-output form  
%
%       yr = Mru*u + Mrd*d + Mrf*f + Mrw*w ,
%
%  with the Laplace- or Z-transformed reference model outputs yr, control  
%  inputs u, disturbance inputs d, fault inputs f, and noise inputs w, and 
%  with Mru, Mrd, Mrf, and Mrw the corresponding transfer-function matrices.  
%  The inputs u, d, f, and w of SYSR correspond to four input groups 
%  named, respectively, {'controls','disturbances','faults','noise'}. 
%  Any of the input groups can be void, in which case, the corresponding 
%  transfer-function matrix is assumed to be zero. 
%
%  The stable filter Q, determined in a standard state-space form, 
%  corresponds to the input-output (implementation) form
%
%            r = Q*[ y ] = Qy*y + Qu*u .
%                  [ u ]
%
%  The inputs y and u of the resulting filter Q are grouped in 
%  two input groups {'outputs','controls'}, respectively, and the output 
%  group {'residuals'} is defined for the residuals r. 
%
%  The filter Q is the solution of the EMMP and is determined to satisfy
%  the linear rational equation
%
%     Q*[ Gu Gd Gf Gw ] = M*[ Mru Mrd Mrf Mrw ] ,                 (1)
%       [ I  0  0  0  ]
%
%  where M is a stable, diagonal, and invertible transfer-function matrix 
%  to be determined such that Q is stable. 
%
%  The stable filter R (also called the internal form of Q), represents the 
%  state-space realization of the updated reference filter, i.e., 
%  R = M*[ Mru Mrd Mrf Mrw ] and has the same input groups as SYSR. 
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
%  The filter Q{i} is the solution of an EMMP and is determined to satisfy
%  the linear rational equation
%
%     Q{i}*[ Gu Gd Gf Gw ] = M{i}*[ Mru_i Mrd_i Mrf_i Mrw_i ] ,         (2)
%          [ I  0  0  0  ]
%
%  where M{i} is a stable, diagonal, and invertible transfer-function matrix 
%  to be determined such that Q{i} is stable. The state realizations of 
%  M{1},...,M{nb} are returned in the cell array M.  
%
%  The stable filter R{i}(also called the internal form of Q{i}), 
%  represents the state-space realization of the updated reference filter, 
%  i.e., R{i} = M{i}*[Mru_i Mrd_i Mrf_i Mrw_i] and has the same input 
%  groups as SYSR{i}. 
%
%  The OPTIONS structure allows to specify various user options, as
%  follows:
%  OPTIONS.tol       - relative tolerance for rank computations
%                      (Default: internally determined value)
%  OPTIONS.tolmin    - absolute tolerance for observability tests
%                      (Default: internally determined value)
%  OPTIONS.smarg     - prescribed stability margin for the poles of the 
%                      filters Q and R
%                      (Default: -sqrt(eps) for a continuous-time system 
%                               and 1-sqrt(eps) for a discrete-time system) 
%  OPTIONS.sdeg      - prescribed stability degree for the poles of the 
%                      filters Q and R
%                      (Default:  -0.05 for a continuous-time system and
%                                  0.95  for a discrete-time system) 
%  OPTIONS.poles     - specifies a complex conjugated set of desired poles 
%                      within the stability domain to be assigned for the 
%                      filters Q and R (Default: []) 
%  OPTIONS.simple    - option to employ a simple proper basis for synthesis 
%                      true  - use a simple basis; the orders of the  
%                              basis vectors are provided in INFO.deg
%                      false - no simple basis computed (default) 
%  OPTIONS.minimal   - specifies the option to perform least order 
%                      synthesis of the filter Q
%                      true  - perform least order synthesis (default)
%                      false - no least order synthesis performed  
%  OPTIONS.regmin    - specifies the regularization option with least order 
%                      left annihilator selection
%                      true  - perform least order selection (default)
%                      false - no least order selection performed  
%  OPTIONS.tcond     - maximum alowed condition number of the employed 
%                      non-orthogonal transformations (Default: 1.e4).
%  OPTIONS.normalize - specifies the normalization option for the diagonal 
%                      elements of the updating filter M, as follows: 
%                      'gain'    - scale with the gains of the 
%                                  zero-pole-gain representation (default)
%                      'dcgain'  - scale with the DC-gains 
%                      'infnorm' - scale with the values of infinity-norm
%  OPTIONS.freq      - complex frequency value to be employed to check the 
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
%                      This option can be used in conjunction with the no 
%                      least order synthesis option OPTIONS.minimal = false. 
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
%                      stable, diagonal and invertible updating matrix M{i}
%                      in (2). 
%  INFO.freq         - employed frequency FREQ to check left invertibility;
%                      FREQ = [] if no frequency-based left invertibility 
%                      check was performed. 
%  INFO.HDesign      - design matrix H employed for the synthesis of 
%                      the fault detection filter Q. 
%                      If SYSR contains nb filters, then INFO.HDesign is a 
%                      an nb x 1 cell array, with INFO.HDesign{i} 
%                      containing the design matrix H_i employed for the 
%                      synthesis of the filter Q{i}; 
%                      H = [] if no design matrix was involved. 

%  Copyright 2017-2019 A. Varga
%  Author:    A. Varga, 07-04-2017.
%  Revisions: A. Varga, 29-07-2017, 28-08-2017, 02-11-2017, 21-02-2018,
%                       06-07-2018, 26-10-2018, 07-07-2019.
%
%  Method: The synthesis Procedures EMM and EMMS from [1] are implemented.
%  Procedure EMM relies on the model-matching synthesis method proposed in 
%  [2], while Procedure EMMS uses the inversion-based method proposed in [3]. 
%  Procedure EMM is generally employed, unless a strong exact fault 
%  detection and isolation problem (strong EFDIP) is solved, in
%  which case Procedure EMMS is used. 
%
%  If SYSR is the given filter, then the strong EFDIP corresponds to the 
%  choice of SYSR with Mru = 0, Mrd = 0, Mrf invertible, and Mrw = 0. 
%  In this case, only the input group {'faults'} must be specified for
%  SYSR. The solution of a fault estimation problem can be targeted by
%  choosing Mrf = I and checking that the resulting INFO.M = I. 
%
%  If SYSR contains a bank of given filters, then the strong EFDIP 
%  corresponds to the choice of SYSR{i} with Mru_i = 0, Mrd_i = 0, 
%  Mrw_i = 0 and Mrf_i set as the i_th row of the mf x mf identity marix.
%  In this case, only the input group {'faults'} must be specified for each
%  SYSR{i}. The solution of a fault estimation problem can be targeted by 
%  choosing Mrf_i as above and checking that the resulting INFO.M{i} = I. 
%
%  References:
%  [1] A. Varga, "Solving Fault Diagnosis Problems - Linear Synthesis 
%      Techniques", Springer Verlag, 2017; sec. 5.6.
%  [2] A. Varga, "New computational approach for the design of fault
%      detection and isolation filters". 
%      In M. Voicu (Ed.), "Advances in Automatic Control", vol. 754 of 
%      The Kluwer International Series in Engineering and Computer Science, 
%      Kluwer Academic Publishers, pp. 367-381, 2003.
%  [3] A. Varga. "New computational paradigms in solving fault detection 
%      and isolation problems". 
%      Annual Reviews in Control, 37:25–42, 2013. 

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
else
    smax = -sqrt(eps);  smin = -inf;
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
   end
   for i = 1:nb
       if Ts ~= sysr{i}.Ts
          error('All components of the reference model SYSR must have the same sampling time')
       end
       if ~isempty(sysr{i}) && ~isstable(sysr{i})
          error('All components of the reference model SYSR must be stable systems')
       end
   end
   rdim(i) = size(sysr{i},1); % number of i-th residual outputs
else
   error('SYSR must be either a cell array or a LTI state space object')
end

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
   normalize = 'gain';
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
   freq = [];
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
          if ~isempty(HDesign{i})
              if size(HDesign{i},1) ~= rdim(i)
                 error(['Row dimension of OPTIONS.HDesign{',num2str(i),'} must be equal to ',num2str(rdim(i))])
              end
              if size(HDesign{i},1) ~= rank(HDesign{i})
                 error(['OPTIONS.HDesign{',num2str(i),'} must have full row rank'])
              end
          end
      end
   end
   emptyHD = all(cellfun('isempty',HDesign));
else
   HDesign = cell(1,nb);
   emptyHD = true;
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

% set options for LCF-based stabilization to be used for final synthesis
opts_glcf = struct('tol',tol,'tolmin',tolmin, 'mininf',true,...
                   'mindeg',true,'sdeg',sdeg,'smarg',smarg,'poles',poles);
% set options for solving linear rational matrix equations               
opts_glsol = struct('tol',1.e-7,'mindeg',minimal);

Q = cell(nb,1); R = cell(nb,1); M = cell(nb,1); Htemp = cell(nb,1);
if all(cellfun('isempty',inpru)) && all(cellfun('isempty',inprd)) && ~minimal
   % perform either Procedure EMM or EMMS if SYSR has no 'controls' and
   % 'disturbances' input groups and no least order option is selected 
      
   % 1) Compute Q1, the left nullspace of [Gu Gd;I 0], and  
   %    R = [Rf1 Rw1] = Q1*[Gf Gw;0 0]
   %    Q_Rfw contains Q1 in Q_Rfw(:,1:p+mu) and R in Q_Rfw(:,p+mu+1:end),
   %    where Rf1 is in Q_Rfw(:,p+mu+1:p+mu+mf) and 
   %    Rw1 is in Q_Rfw(:,p+mu+mf+1:end)

   % set options for nullspace computation
   opts_glnull = struct('tol',tol,'m2',mf+mw,'simple',simple);
   % form [ Gu Gd Gf Gw; I 0 0 0 ] 
   syse = [sysf(:,[inpu inpd inpf inpw]); eye(mu,m)];
   [QR,info1] = glnull(syse,opts_glnull); 
   
   nvec = size(QR,1);       % number of basis vectors
   % check solvability conditions
   if nvec == 0,
      error('Empty nullspace basis: the EBMMP is not solvable')
   end
   degs = info1.degs;       % degrees of a minimal polynomial basis
   tcond1 = info1.tcond;    % condition number of employed transformations
   
   indf = p+mu+(1:mf);      % input indices of Rf in QR
   indfw = p+mu+(1:mf+mw);  % input indices of [ Rf Rw ] in QR
   Rfw = QR(:,indfw);

   
   if isempty(freq)
      freq = rand; 
   end
   Rftest = evalfr(QR(:,indf),freq);
   for ib = 1:nb
     rdimi = rdim(ib);
     if mf+mw > 0
       % set H for checking the solvability condition
       if emptyHD
          Htemp{ib} = eye(nvec);
       else
          degs = [];
          [mh,nh] = size(HDesign{ib});
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
       % this is the usual case with nonempty [Rf1 Rw1]
       if rdimi == mf && nvec >= mf && ...
         ((tol > 0 && rank(evalfr(sysr{i}(:,'faults'),freq),tol) == mf) || ...
          (tol == 0 && rank(evalfr(sysr{i}(:,'faults'),freq)) == mf)) && ...
         ((tol > 0 && rank(Htemp{ib}*Rftest,tol) == mf) || ...
          (tol == 0 && rank(Htemp{ib}*Rftest) == mf)) 
         % perform Procedure EMMS if Gw = 0, Rf is invertible,
         % and the reduced Rf1 is left invertible
         % flip degrees of a minimal polynomial basis
         if ~simple, degs = flip(degs); end
         finish = (nvec == mf);  % set termination flag
         nq = order(QR);
         if ~finish
            nout = mf;       % initialize number of selected basis vectors
            if ~simple && regmin
               % permute states to speedup glmcover1 
               QR = xperm(QR,nq:-1:1);  
            end
         else
            if emptyHD 
               h = eye(nvec);
            end
         end
         while ~finish     
            % choose nout basis vectors, which potentially lead to a 
            % least order filter with mf outputs:
            % basesel(i,:) contains the indices of candidate basis vectors;
            % ordsel(i)    contains the presumably achievable least orders
            [basesel,ordsel] = emmbasesel(Htemp{ib}*Rftest,degs,nout,simple,tol); 
            % update the synthesis using the selections of candidate 
            % vector(s), starting with the least (potentially) achievable order
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
                               warning('EMMSYN: Expected reduced order not achieved')
                            end
                            h = Htemp{ib}(ip(1:rdimi),:);
                         else
                             [QRfwtest,info2] = glmcover1([Htemp{ib}; eye(nvec)]*QR(ip,:),rdimi,tol);
                         end
                      else   
                         % the case rdim < nout can only happen if no
                         % HDesign is explicitly provided
                         Htemp{ib} = blkdiag(hbase,eye(nvec-nout)); 
                         [QRfwtest,info2] = glmcover1([Htemp{ib}; eye(nvec)]*QR(ip,:),rdimi,tol);
                      end
                   else
                      % here only the case rdim = nout can happen
                      if emptyHD
                         h = Htemp{ib}(ip(1:rdimi),:);
                         QRfwtest = gir(QR(baseind,:),tol,'finite');
                      else
                         QRfwtest = gir(Htemp{ib}*QR,tol,'finite');
                      end
                   end
                end
                % check invertibility of compressed Rf1; 
                if ~simple && regmin
                   % dismiss minimal design if the check fails
                   Rtest1 = evalfr(QRfwtest(:,indf),freq); 
                   if (tol > 0 && rank(Rtest1,tol) == mf) || ...
                      (tol == 0 && rank(Rtest1) == mf)
                      % adjust condition number of employed transformations
                      tcond1 = max([tcond1; info2.fnorm;info2.tcond]);
                      if tcond1 > tcond
                         disp(['EMMSYN: Possible loss of numerical stability',...
                               ' due to ill-conditioned transformations'])
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
            if nout > nvec && ~finish
               error('Something wrong: try perheps with another test frequency')
            end
         end
         
         if emptyHD
            Htemp{ib} = h;
         end
         % compute the irreducible realization of Qtilde = Rf*(inv(Rf1)*Q1)  
         % by first solving the linear rational matrix equation Rf*Q2 = Q1
         Q2 = grsol(QR(:,[indf 1:p+mu]),p+mu,struct('tol',tol));
         Qtilde = gminreal(sysr{ib}(:,'faults')*Q2,tol);       
       else
         Htemp{ib} = []; 
         freq = [];
         if isempty(inprw{ib}) 
            [Q2,info2] = glsol(Rfw(:,1:mf),sysr{ib},opts_glsol);
         else
            [Q2,info2] = glsol(Rfw,sysr{ib},opts_glsol);
         end
%            
%          % perform Procedure EMM in the general case
%          % 2) Solve Q2*Rfw = [Rf Rw] and form Q = Q2*Q1 .
%          Rref = sysr;
%          if isempty(inprf) && ~isempty(inpf)
%             Rref = [zeros(rdim,mf) Rref ];
%          end
%          if isempty(inprw) && ~isempty(inpw)
%             Rref = [Rref zeros(rdim,mw)];
%          end        
%          [Q2,info2] = glsol(Rfw,Rref,opts_glsol);
         Qtilde = Q2*QR(:,1:p+mu);
         tcond1 = max([tcond1; info2.fnorm;info2.tcond]);
       end
        
       % 3) compute diagonal M such that M*Q has a desired stability degree;  
       %    update Q <- M*Q
       Q{ib} = ss(zeros(rdimi,p+mu)); M{ib} = ss(zeros(rdimi));
       for i=1:rdimi
          % LCF is performed after removing unobservable eigenvalues
          [Qi,Mi] = glcf(gir(Qtilde(i,:),tolmin,'obs'),opts_glcf);
          switch job
             case 0  % gain
                sc = get(zpk(Mi),'k');  % scale with gain 
             case 1  % dcgain
                sc = dcgain(Mi);        % scale with dcgain 
             case 2  % infnorm
                sc = norm(Mi,inf);      % scale with infinity norm
          end          
          %Q = [Q;Qi/sc]; M = append(M,Mi/sc);
          Q{ib}(i,:) = Qi/sc; M{ib}(i,i) = Mi/sc;
       end
     else
      % this is the case with [Rf1 Rw1] empty; M must not be diagonal
      % 3) compute M such that M*Q has a desired stability degree;  
      %    update Q <- M*Q
      [Q{ib},M{ib}] = glcf(QR,opts_glcf); 
     end
     R{ib} = M{ib}*sysr{ib};     % R = M*SYSR    
     if ~isempty(inpw) && isempty(inprw{ib})
        R{ib} = gir([R{ib} Q{ib}*[sysf(:,inpw); zeros(mu,mw)]],tol); 
     end
   end
else
   % apply the two-step procedure to solve the EMMP if the least order
   % synthesis option has been selected or SYSR has either 
   % 'controls' or 'disturbances' input groups, or both
   degs = []; 
   Htemp = cell(nb,1);
   for ib = 1:nb
     rdimi = rdim(ib);
     if isempty(inprw{ib})
        % apply the two-step procedure for the case
        % form Ge = [ Gu Gd Gf; I 0 0 ] 
        m1 = m-mw;
        syse = [sysf(:,[inpu inpd inpf]); eye(mu,m1)];
        rinp = zeros(0,m1);
        if ~isempty(inpru{ib}) 
           rinp = [rinp;eye(mu,m1)];
        end
        if ~isempty(inprd{ib}) 
           rinp = [rinp; zeros(md,mu) eye(md,m1-mu)];
        end
        if ~isempty(inprf{ib}) 
           rinp = [rinp; zeros(mf,mu+md) eye(mf,m1-mu-md)];
        end
        % form explicitly Rref = [ Mru_i Mrd_i Mrf_i ] 
        Rref = sysr{ib}*rinp;
     else
        % apply the two-step procedure for the general case
        % form Ge = [ Gu Gd Gf Gw; I 0 0 0 ] 
        syse = [sysf(:,[inpu inpd inpf inpw]); eye(mu,m)];
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
        % form explicitly Rref = [Mru_i Mrd_i Mrf_i Mrw_i ] 
        Rref = sysr{ib}*rinp;
     end
   
     % 1) Solve Q*Ge = Rref .
     [Qtilde,info2] = glsol(syse,Rref,opts_glsol);
     tcond1 = max([info2.fnorm;info2.tcond]);
   
     % 2) compute diagonal M such that M*Q has a desired stability degree;  
     %    update Q <- M*Q
     Q{ib} = ss(zeros(rdimi,p+mu)); M{ib} = ss(zeros(rdimi));
     for i=1:rdimi
       % LCF is performed after removing unobservable eigenvalues
       [Qi,Mi] = glcf(gir(Qtilde(i,:),tolmin,'obs'),opts_glcf);
       switch job
          case 0  % gain
             sc = get(zpk(Mi),'k');  % scale with gain 
          case 1  % dcgain
             sc = max(1.e-4,dcgain(Mi));        % scale with dcgain 
          case 2  % infnorm
             sc = norm(Mi,inf);      % scale with infinity norm
       end          
       %Q = [Q;Qi/sc]; M = append(M,Mi/sc);
       Q{ib}(i,:) = Qi/sc; M{ib}(i,i) = Mi/sc;
     end
     R{ib} = M{ib}*sysr{ib}; 
     if ~isempty(inpw) && isempty(inprw{ib})
        R{ib} = gir([R{ib} Q{ib}*[sysf(:,inpw); zeros(mu,mw)]],tol); 
     end
   end
end


for ib = 1:nb
  rdimi = rdim(ib);
  % transform to standard state-space
  Qib = gss2ss(Q{ib},tol);

  % set output variables
  set(Qib,'InputGroup',struct('outputs',1:p,'controls',p+(1:mu)));
  set(Qib,'OutputGroup',struct('residuals',1:rdimi));
  Q{ib} = Qib;

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
     if ~isempty(inprf{ib})
        Rib.InputGroup.faults = ioff+(1:mf);
        ioff = ioff+mf;
     end
     if ~isempty(inprw{ib})
        Rib.InputGroup.noise = ioff+mf+(1:mw);
     end
     R{ib} = Rib;
  end
end

if sysrss
   Q = Q{1}; R = R{1}; 
end

if nargout > 2
   if sysrss
       info = struct('tcond',tcond1,'degs',degs,'M',M{1},'freq',freq,'HDesign',Htemp{1}); 
   else
       info = struct('tcond',tcond1,'degs',degs,'M',{M},'freq',freq,'HDesign',{Htemp}); 
   end
end

% end EMMSYN
end

