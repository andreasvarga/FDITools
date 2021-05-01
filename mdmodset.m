function sysm = mdmodset(sys,inputs)
%MDMODSET Setup of models for solving model detection synthesis problems
%  SYSM = MDMODSET(SYS,INPUTS) builds for a multiple model SYS with N 
%  component models, where the j-th system has the input-output form 
%
%            y_j = G_j*u0_j
%
%  with the output vector y_j and input vector u0_j, a multiple model SYSM 
%  with N component models, with the j-th model in the corresponding 
%  input-output form
%
%      y_j = G_{u,j}*u + G_{d,j}*d_j + G_{w,j}*w_j ,
%
%  where the control inputs u, disturbance inputs d_j and noise inputs w_j
%  are derived from the system input vector u0_j in accordance 
%  with the input groups specified in the INPUTS structure. 
%  The variables y_j, u0_j, u, d_j, and w_j are Laplace- or Z-transformed
%  vectors, and G_j, G_{u,j}, G_{d,j} and G_{w,j} are the corresponding 
%  transfer-function matrices. 
%  The inputs u, d_j and w_j of the j-th model of SYSM correspond to three
%  input groups named, respectively, {'controls','disturbances','noise'}.
%  Any of these input groups may be void. 
%  The resulting multiple model SYSM preserves the common sampling time 
%  of the original component models of SYS. 
%
%  SYS contains the multiple model specified either as an N-vector of 
%  LTI systems or a cell array with N elements, whose components are 
%  LTI systems with the same number of outputs and control inputs, and 
%  the same sampling time. Correspondingly, the resulting SYSM is an
%  N-vector of LTI systems or a cell array of N LTI systems, respectively. 
%
%  Note: Excepting the sampling time, no other model attributes are
%  inherited by the component models of SYSM from the original component 
%  models of SYS. 
%  
%  The INPUTS structure allows to define various input groups, as
%  follows:
%  INPUTS.controls     - indices of the control inputs (Default: void)
%  INPUTS.c            - alternative short form to specify the indices of 
%                        the control inputs (Default: void)
%  INPUTS.disturbances - indices of the disturbance inputs or an 
%                        N-dimensional cell array, with 
%                        INPUTS.disturbances{j} containing the indices of 
%                        the disturbance inputs of the j-th component 
%                        system (Default: void)
%  INPUTS.d            - alternative short form to specify the indices of  
%                        the disturbance inputs or an N-dimensional cell 
%                        array, with INPUTS.d{j} containing the indices of 
%                        the disturbance inputs of the j-th component 
%                        system (Default: void)
%  INPUTS.noise        - indices of the noise inputs or an 
%                        N-dimensional cell array, with 
%                        INPUTS.noise{j} containing the indices of 
%                        the noise inputs of the j-th component 
%                        system (Default: void)
%  INPUTS.n            - alternative short form to specify the indices of  
%                        the noise inputs or an N-dimensional cell 
%                        array, with INPUTS.n{j} containing the indices of 
%                        the noise inputs of the j-th component 
%                        system (Default: void)

%  Copyright 2018 A. Varga
%  Author:    A. Varga, 01-06-2018.
%  Revisions:  


narginchk(1,2)
nargoutchk(0,1)

% check input system form
if isa(sys,'cell')
   validateattributes(sys, {'cell'},{'vector'},'','SYS')
   N = length(sys);
   for i = 1:N
       if ~isa(sys{i},'ss')
          error('The input SYS must be a cell array of LTI state space objects')
       end
   end
   Ts = sys{1}.Ts;
   m = zeros(N,1);
   % determine systems dimensions
   [p,m(1)] = size(sys{1});
   for i = 2:N
       [noi,m(i)] = size(sys{i}); 
       if noi ~= p 
          error('All component models must have the same number of outputs')
       end
       if Ts ~= sys{i}.Ts
          error('All component models must have the same sampling time')
       end
   end
elseif isa(sys,'ss')
   % determine system dimensions
   [p,m,M,N] = size(sys);
   if min(M,N) ~= 1
       error('Only one-dimensional arrays of systems are supported')
   end
   Ts = sys.Ts; 
else
   error('SYS must be either a cell array or a LTI state space object')
end

mmin = min(m);

if nargin > 1
   validateattributes(inputs,{'struct'},{'nonempty'},'','INPUTS')
else
   % no input groups defined
   if isa(sys,'cell')
       sysm = cell(N,1);
       for j = 1:N
           sysm{j} = ss(zeros(p,0));
           sysm{j}.Ts = Ts;
       end
   else
       sysm = ss(zeros(p,0,M,N));
       sysm.Ts = Ts;
   end
   return
end
    
% decode input information
% control inputs
if isfield(inputs,'controls')
    inpu = inputs.controls;
    if ~isempty(inpu)
       validateattributes(inpu, {'double'},{'integer','vector','>',0,'<=',mmin},'','INPUTS.controls') 
    end
    mu = length(inpu);  
else
    inpu = []; mu = 0;
end
if isempty(inpu) 
    if isfield(inputs,'c')
       inpu = inputs.c;
       if ~isempty(inpu)
          validateattributes(inpu, {'double'},{'integer','vector','>',0,'<=',mmin},'','INPUTS.c') 
       end
       mu = length(inpu);  
    end
else
    if isfield(inputs,'c')
       error('Either INPUTS.controls or INPUTS.c can be only specified')
    end
end

% disturbance inputs
if isfield(inputs,'disturbances') 
   inpd = inputs.disturbances;
else
   inpd = [];
end
if isfield(inputs,'d')
   if isempty(inpd)  
      inpd = inputs.d;
   else
      error('Either INPUTS.disturbances or INPUTS.d can be only specified')
   end
end
if ~isempty(inpd)
    if isa(inpd,'cell');
       Nd = length(inpd);
       validateattributes(inpd, {'cell'},{'vector'},'','INPUTS.disturbances or INPUTS.d')
       if N ~= Nd
          error('The number of index sets in INPUTS.disturbances or INPUTS.d must be equal to the number of systems')
       end
       md = zeros(N,1);
       for j = 1:N
           validateattributes(inpd{j}, {'double'},{'integer','vector','>',0,'<=',mmin},'','INPUTS.disturbances or INPUTS.d')
           md(j) = length(inpd{j});
       end
    else
       if ~isempty(inpd)
           validateattributes(inpd, {'double'},{'integer','vector','>',0,'<=',mmin},'','INPUTS.disturbances or INPUTS.d') 
       end
       if isa(sys,'cell')
          inpdsav = inpd;
          inpd = cell(N,1);
          for j = 1:N, inpd{j} = inpdsav; end
          md = length(inpdsav)*ones(N,1);
       else
          md = length(inpd);
       end
    end
else
    if isa(sys,'cell')
       inpd = cell(N,1); md = zeros(N,1);
    else
       inpd = []; md = 0;
    end
end

% noise input
if isfield(inputs,'noise') 
   inpw = inputs.noise;
else
   inpw = [];
end
if isfield(inputs,'n')
   if isempty(inpw)  
      inpw = inputs.n;
   else
      error('Either INPUTS.noise or INPUTS.n can be only specified')
   end
end
if ~isempty(inpw)
    if isa(inpw,'cell') 
       Nw = length(inpw);
       validateattributes(inpw, {'cell'},{'vector'},'','INPUTS.noise or INPUTS.n')
       if N ~= Nw
          error('The number of index sets in INPUTS.noise or INPUTS.n must be equal to the number of systems')
       end
       mw = zeros(N,1);
       for j = 1:N
           validateattributes(inpw{j}, {'double'},{'integer','vector','>',0,'<=',mmin},'','INPUTS.noise or INPUTS.n')
           mw(j) = length(inpw{j});
       end
    else
       if ~isempty(inpw)
          validateattributes(inpw, {'double'},{'integer','vector','>',0,'<=',mmin},'','INPUTS.noise or INPUTS.n') 
       end
       if isa(sys,'cell')
          inpwsav = inpw;
          inpw = cell(N,1);
          for j = 1:N, inpw{j} = inpwsav; end
          mw = length(inpwsav)*ones(N,1); 
       else
          mw = length(inpw); 
       end
   end
else
    if isa(sys,'cell')
       inpw = cell(N,1); mw = zeros(N,1);
    else
       inpw = []; mw = 0;
    end
end

% setup of the synthesis multiple model
if isa(sys,'ss')
   sysm = ss(zeros(p,mu+md+mw,M,N));
   Ts = sys.Ts;
   sysm.Ts = Ts;
   for im = 1:M
       for jm = 1:N
           [A,B,C,D,E,Ts] = dssdata(sys,im,jm);
            n = size(A,1);
            standsys = isequal(E,eye(n));
            Be = [B(:,inpu) B(:,inpd) B(:,inpw)];
            De = [D(:,inpu) D(:,inpd) D(:,inpw)];
            if standsys
               sysm(:,:,im,jm) = ss(A,Be,C,De,Ts);
            else
               sysm(:,:,im,jm) = dss(A,Be,C,De,E,Ts);
            end
       end
   end
   % set input groups
   InpG = struct('controls',1:mu,...
              'disturbances',mu+(1:md),...
              'noise',mu+md+(1:mw));
   set(sysm,'InputGroup',InpG);
else
   sysm = cell(N,1);
   for jm = 1:N
       [A,B,C,D,E,Ts] = dssdata(sys{jm});
       n = size(A,1);
       standsys = isequal(E,eye(n));
       Be = [B(:,inpu) B(:,inpd{j}) B(:,inpw{j})];
       De = [D(:,inpu) D(:,inpd{j}) D(:,inpw{j})];
       if standsys
          syst = ss(A,Be,C,De,Ts);
       else
          syst = dss(A,Be,C,De,E,Ts);
       end
       % set input groups
       InpG = struct('controls',1:mu,...
                     'disturbances',mu+(1:md(jm)),...
                     'noise',mu+md(jm)+(1:mw(jm)));
       set(syst,'InputGroup',InpG);
       sysm{jm} = syst;
   end
end
% end MDMODSET
end
