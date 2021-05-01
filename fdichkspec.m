function [rdims,orders,leastorders] = fdichkspec(sysf,SFDI,options) 
%FDICHKSPEC  Feasibility analysis of a set of FDI specifications.
%  [RDIMS,ORDERS,LEASTORDERS] = FDIGENSPEC(SYSF,SFDI,OPTION) determines for 
%  a LTI system SYSF with mf fault inputs and a given logical N x mf matrix 
%  SFDI containing N FDI specifications, the N-dimensional vector RDIMS,
%  whose i-th nonzero component RDIMS(i), contains the number of residual 
%  outputs of a minimal nullspace basis based FDI filter which can be used 
%  to achieve the i-th specification contained in SFDI(i,:). If SFDI = []
%  or not specified, then the fault inputs are ignored.  
%  If RDIMS(i) = 0, the i-th specification is not feasible. 
%  ORDERS is an N-dimensional vector, whose i-th component ORDERS(i) 
%  contains for a feasible specification SFDI(i,:) the order of the minimal 
%  nullspace basis based FDI filter (see above). If the i-th specification
%  is not feasible, then ORDERS(i) is set to -1.
%  LEASTORDERS is an N-dimensional vector, whose i-th component 
%  LEASTORDERS(i) contains, for a feasible specification SFDI(i,:), the 
%  least achievable order for a scalar output FDI filter which can be used  
%  to achieve the i-th specification. If the i-th specification is not 
%  feasible, then LEASTORDER(i) is set to -1.
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
%  The OPTIONS structure allows to specify various user options, 
%  as follows:
%  OPTIONS.tol       - tolerance for rank determinations
%                      (Default: internally computed)
%  OPTIONS.tolmin    - absolute tolerance for observability tests
%                      (Default: internally determined value)
%  OPTIONS.FDTol     - threshold for assessing weak specifications
%                      (see function FDITSPEC) (Default: 0.0001)
%  OPTIONS.FDGainTol - threshold for assessing strong specifications,
%                      i.e., threshold for nonzero frequency-responce
%                      gains for all frequency values specified in
%                      OPTIONS.FDFreq (see function FDISSPEC) 
%                      (Default: 0.01)
%  OPTIONS.FDFreq    - vector of real frequency values for strong  
%                      detectability checks (Default: [])
%
%  See also FDITSPEC and FDISSPEC.

%  Copyright 2018 A. Varga
%  Author:      A. Varga, 24-10-2018.
%  Revision(s): 
%
%  Method: The nullspace method of [1] is successively employed to 
%  determine FDI filters as minimal left nullspace bases which solve 
%  suitably formulated fault detection problems. 
%
%  References:
%  [1] Varga A.
%      On computing nullspace bases – a fault detection perspective. 
%      Proc. IFAC 2008 World Congress, Seoul, Korea, pages 6295–6300, 2008.

narginchk(1,3)
nargoutchk(0,3)


% check input system form
if ~isa(sysf,'ss')
   error('The input system SYSF must be an SS object')
end

if nargin < 2
   SFDI = [];
end

if nargin < 3
   options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

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

if isempty(SFDI)
   SFDI = false(0,mf);
   N = 1;
else
   validateattributes(SFDI, {'logical'},{'binary'},'','OPTIONS.SFDI') 
   if size(SFDI,2) ~= mf
      error('Structure matrix incompatible with the number of faults')
   end
   N = size(SFDI,1);      % number of specifications 
end


m = mu+md+mf;
p = size(sysf,1);
lfreq = length(FDFreq);     % number of frequency values
rdims = zeros(N,1);
orders = -ones(N,1);
leastorders = -ones(N,1);


% set options for nullspace computation
if strongFD
   % set default stability degree
   if sysf.Ts ~= 0
      sdegdefault = 0.95;
   else
      sdegdefault = -0.05;
   end
   % set options for nullspace computation with stabilization
   opts_glnull = struct('tol',tol,'sdeg',sdegdefault);
else 
   opts_glnull = struct('tol',tol);
end
% set options for LCF-based stabilization to be used for solvability checks
opts_glcf_default = struct('tol',tol,'tolmin',tolmin);
for i = 1:N
   if isempty(SFDI)
      syse = [sysf(:,[inpu inpd]); eye(mu,mu+md)];
      opts_glnull.m2 = 0;
      indff = [];
   else
      % set input groups for first design
      inddf = (SFDI(i,:) == 0); indff = find(SFDI(i,:) ~= 0);
      syse = [sysf(:,[inpu inpd inpf(inddf) inpf(~inddf)]); eye(mu,m)];
      m2f = length(indff);
      opts_glnull.m2 = m2f;
   end
   %
   % compute a left nullspace basis Q of G = [Gu Gd Gdf; I 0 0] = 0 
   % obtain QR = [ Q Rf ], where Rf = Q*[Gf;0] 
   [QR,info] = glnull(syse,opts_glnull); 
   nvec = size(QR,1);          % number of basis vectors
   % check solvability conditions
   if nvec == 0 
      continue
   end
   
   nq = order(QR);            % order of minimal basis
   degs = flip(info.degs);    % degrees of a minimal polynomial basis
   if isempty(indff)
      rdims(i) = nvec;
      orders(i) = nq;
      leastorders(i) = min(degs);
      continue
   end
   
   indf = p+mu+(1:m2f);
   feasible = true;
   if strongFD 
      S = fdisspec(QR(:,indf),FDGainTol,FDFreq);
      % check strong detectability conditions 
      for ii = 1:lfreq
          if ~all(max(S(:,:,ii),[],1))
             feasible = false;
             break
          end  
      end
   else
      % check weak detectability conditions 
      S = fditspec(QR(:,indf),tol,FDTol);
      if ~all(max(S,[],1))
         feasible = false;
      end   
   end
   if feasible
      rdims(i) = nvec;
      orders(i) = nq; 
      leastorders(i) = nq;
      if nargout > 2
         if nvec > 1
            finish = false;    % set termination flag
            nout = 1;          % initialize number of selected basis vectors
            QR = xperm(QR,nq:-1:1);  % permute states to speedup glmcover1
            while ~finish     
               % choose nout basis vectors, which potentially lead to a least order
               % filter with rdim outputs:
               % basesel(i,:) contains the indices of candidate basis vectors;
               % ordsel(i)    contains the presumably achievable least orders
               [basesel,ordsel] = efdbasesel(S,degs,1,nout,false); 
               %
               % update the synthesis using the selections of candidate vector(s),
               % starting with the least (potentially) achievable order
               for ibas = 1:size(basesel,1);
                  baseind = basesel(ibas,:); % indices of current basis selection
                  ip = [baseind, setdiff(1:nvec,baseind)]; 
                  if nout == 1
                     QRfwtest = glmcover1(QR(ip,:),1,tol);
                     if order(QRfwtest) ~= ordsel(ibas)
                        warning('FDICHKSPEC: Expected reduced order not achieved')
                     end
                  else  
                     % build a linear combination of the first nout vectors 
                     h = [ rand(1,nout) zeros(1,nvec-nout) ]; 
                     QRfwtest = glmcover1([h; eye(nvec)]*QR(ip,:),1,tol);
                  end
                  % check complete fault detectability of the current design; 
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
                     finish = true;
                     leastorders(i) = order(QRfwtest);
                     break
                  end
               end
               nout = nout+1;
               if nout > nvec
                 finish = true;
               end
            end
         end
      end
   end
end

% end FDICHKSPEC  
end
