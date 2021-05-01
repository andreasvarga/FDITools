function sysf = fdimodset(sys,inputs)
%FDIMODSET Setup of models for solving FDI synthesis problems
%  SYSF = FDIMODSET(SYS,INPUTS) builds for the LTI state-space system 
%  SYS with an input-output form 
%
%            y = G*u0
%
%  with the output vector y and input vector u0, the LTI state-space system
%  SYSF with additive faults with the corresponding input-output form
%
%            y = Gu*u + Gd*d + Gf*f + Gw*w + Ga*aux ,
%
%  where the control input vector u, disturbance input vector d, 
%  fault input vector f, noise input vector w, and auxiliary input 
%  vector aux, are derived from the system input vector u0 in accordance 
%  with the input groups specified in the INPUTS structure. 
%  The variables y, u0, u, d, f, w and aux are Laplace- or Z-transformed
%  vectors, and G, Gu, Gd, Gf, Gw, and Ga are the corresponding 
%  transfer-function matrices. 
%  The inputs u, d, f, w and aux of SYSF correspond to five input groups 
%  named, respectively, {'controls','disturbances','faults','noise','aux'}.
%  Any of these input groups may be void. 
%  The resulting model SYSF inherits the sampling time of the 
%  original model SYS. 
%
%  FDIMODSET can also be employed if SYS is an N x 1 or 1 x N array of LTI 
%  models, in which case the resulting SYSF is also an N x 1 or 1 x N array
%  of LTI models, respectively.
%
%  Note: Excepting the sampling time, no other model attributes are
%  inherited by SYSF from the original model SYS. 
%  
%  The INPUTS structure allows to define various input groups, as
%  follows:
%  INPUTS.controls     - indices of the control inputs (Default: void)
%  INPUTS.c            - alternative short form to specify the indices of 
%                        the control inputs (Default: void)
%  INPUTS.disturbances - indices of the disturbance inputs (Default: void)
%  INPUTS.d            - alternative short form to specify the indices of 
%                        the disturbance inputs (Default: void)
%  INPUTS.faults       - indices of the fault inputs (Default: void)
%  INPUTS.f            - alternative short form to specify the indices of 
%                        the fault inputs (Default: void)
%  INPUTS.faults_sen   - indices of the outputs subject to sensor faults
%                        (Default: void)
%  INPUTS.fs           - alternative short form to specify the indices of  
%                        the outputs subject to sensor faults 
%                        (Default: void)
%  INPUTS.noise        - indices of the noise inputs (Default: void)
%  INPUTS.n            - alternative short form to specify the indices of 
%                        the noise inputs (Default: void)
%  INPUTS.aux          - indices of the auxiliary inputs (Default: void)

%  Copyright 2018 A. Varga
%  Author:    A. Varga, 04-06-2018.
%  Revisions:  


narginchk(1,2)
nargoutchk(0,1)

% check input system form
if ~isa(sys,'ss')
   error('The input system SYS must be an SS object')
end

% determine system dimensions
[p,m,M,N] = size(sys); 
if min(M,N) ~= 1
    error('Only one-dimensional arrays of systems are supported')
end

if nargin > 1
   validateattributes(inputs,{'struct'},{'nonempty'},'','INPUTS')
else
   % no input groups defined
   sysf = ss(zeros(p,0,M,N));
   sysf.Ts = sys.Ts;
   return
end
    



% decode input information
% control inputs
if isfield(inputs,'controls')
    inpu = inputs.controls;
    if ~isempty(inpu)
       validateattributes(inpu, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.controls') 
    end
    mu = length(inpu);  
else
    inpu = []; mu = 0;
end
if isempty(inpu) 
    if isfield(inputs,'c')
       inpu = inputs.c;
       if ~isempty(inpu)
          validateattributes(inpu, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.c') 
       end
       mu = length(inpu);  
    end
else
    % consistency check 
    if isfield(inputs,'c')
       temp = inputs.c;
       if ~isempty(temp)
          validateattributes(temp, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.c')
       end
       if ~isequal(inpu,temp)
           error('INPUTS.controls and INPUTS.c must have the same content')
       end
    end
end

% disturbance inputs
if isfield(inputs,'disturbances')
    inpd = inputs.disturbances;
    validateattributes(inpd, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.disturbances') 
    md = length(inpd);  
else
    inpd = []; md = 0;
end
if isempty(inpd) 
    if isfield(inputs,'d')
       inpd = inputs.d;
       if ~isempty(inpd)
          validateattributes(inpd, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.d')
       end
       md = length(inpd);  
    end
else
    % consistency check 
    if isfield(inputs,'d')
       temp = inputs.d;
       if ~isempty(temp)
          validateattributes(temp, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.d')
       end
       if ~isequal(inpd,temp)
           error('INPUTS.disturbances and INPUTS.d must have the same content')
       end
    end
end

% fault inputs
if isfield(inputs,'faults')
   inpf = inputs.faults;
   if ~isempty(inpf)
      validateattributes(inpf, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.faults') 
   end
   mf = length(inpf);  
else
    inpf = []; mf = 0;
end
if isempty(inpf) 
   if isfield(inputs,'f')
      inpf = inputs.f;
      if ~isempty(inpf)
         validateattributes(inpf, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.f') 
      end
      mf = length(inpf);  
   end
else
    % consistency check 
    if isfield(inputs,'f')
       temp = inputs.f;
       if ~isempty(temp)
          validateattributes(temp, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.f')
       end
       if ~isequal(inpf,temp)
           error('INPUTS.faults and INPUTS.f must have the same content')
       end
    end
end
% noise input
if isfield(inputs,'noise')
   inpw = inputs.noise;
   validateattributes(inpw, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.noise') 
   mw = length(inpw);  
else
   inpw = []; mw = 0;
end
if isempty(inpw) 
   if isfield(inputs,'n')
      inpw = inputs.n;
      if ~isempty(inpw)
         validateattributes(inpw, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.n') 
      end
      mw = length(inpw);  
   end
else
    % consistency check 
    if isfield(inputs,'n')
       temp = inputs.n;
       if ~isempty(temp)
          validateattributes(temp, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.n')
       end
       if ~isequal(inpfw,temp)
           error('INPUTS.noise and INPUTS.n must have the same content')
       end
    end
end

% auxiliary inputs
if isfield(inputs,'aux')
    inpaux = inputs.aux;
    if ~isempty(inpaux)
       validateattributes(inpaux, {'double'},{'integer','vector','>',0,'<=',m},'','INPUTS.aux') 
    end
    maux = length(inpaux);  
else
    inpaux = []; maux = 0;
end

% sensor fault inputs
if isfield(inputs,'faults_sen')
   inpfs = inputs.faults_sen;
   if ~isempty(inpfs)
      validateattributes(inpfs, {'double'},{'integer','vector','>',0,'<=',p},'','INPUTS.faults_sen')
   end
   mfs = length(inpfs);  
else
   inpfs = []; mfs = 0;
end
if isempty(inpfs) 
   if isfield(inputs,'fs')
      inpfs = inputs.fs;
      if ~isempty(inpfs)
         validateattributes(inpfs, {'double'},{'integer','vector','>',0,'<=',p},'','INPUTS.fs')
      end
      mfs = length(inpfs);  
   end
else
    % consistency check 
    if isfield(inputs,'fs')
       temp = inputs.fs;
       if ~isempty(temp)
          validateattributes(temp, {'double'},{'integer','vector','>',0,'<=',p},'','INPUTS.fs')
       end
       if ~isequal(inpfs,temp)
           error('INPUTS.faults_sen and INPUTS.fs must have the same content')
       end
    end
end
if mfs
   Dsf = eye(p); Dsf = Dsf(:,inpfs);
else
   Dsf = zeros(p,0);
end

% setup of the synthesis model
mfe = mf + mfs;
sysf = ss(zeros(p,mu+md+mfe+mw+maux,M,N));
Ts = sys.Ts;
sysf.Ts = Ts;
for im = 1:M
    for jm = 1:N
        [A,B,C,D,E,Ts] = dssdata(sys,im,jm);
         n = size(A,1);
         standsys = isequal(E,eye(n));
         Be = [B(:,inpu) B(:,inpd) B(:,inpf) ...
               zeros(n,mfs) B(:,inpw) B(:,inpaux)];
         De = [D(:,inpu) D(:,inpd) D(:,inpf) ...
               Dsf          D(:,inpw) D(:,inpaux)];
         if standsys
            sysf(:,:,im,jm) = ss(A,Be,C,De,Ts);
         else
            sysf(:,:,im,jm) = dss(A,Be,C,De,E,Ts);
         end
    end
end

% set input groups
InpG = struct('controls',1:mu,...
              'disturbances',mu+(1:md),...
              'faults',mu+md+(1:mfe),...
              'noise',mu+md+mfe+(1:mw),...
              'aux',mu+md+mfe+mw+(1:maux));
set(sysf,'InputGroup',InpG);

% end FDIMODSET
end
