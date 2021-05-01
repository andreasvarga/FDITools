function [dist,fpeak,mind]  = mddist2c(sysm,sys,options) 
% MDDIST2C Computation of distances to a set of component models 
%  [DIST,FPEAK,MIND] = MDDIST2C(SYSM,SYS,OPTIONS) evaluates for a 
%  multiple model SYSM with N component models and a given system model SYS, 
%  the distances DIST between the N component models and SYS.
%  SYSM can be specified either as an N-vector of LTI systems or a 
%  cell array with N elements, whose components are LTI systems with the 
%  same number of outputs and control inputs, and the same sampling time. 
%  The j-th component system is either a continuous- or discrete-time 
%  system in a standard or descriptor state-space form, which corresponds 
%  to the input-output form  
%
%      y_j = Gu_j*u + Gd_j*d_j+ Gw_j*w_j ,
%
%  with the Laplace- or Z-transformed plant outputs y_j, control inputs u, 
%  disturbance inputs d_j and noise inputs w_j, and with Gu_j, Gd_j and 
%  Gw_j the corresponding transfer-function matrices.
%  The control inputs u, common to all component models, must correspond  
%  to the input group named 'controls'. If OPTIONS.cdinp = true, then
%  the same disturbance input d is assumed for all component models 
%  (d_j = d) and must correspond to the input group named 'disturbances'. 
%  SYS is LTI system in state-space form, which corresponds to the 
%  input-output form  
%
%      y = Gu*u + Gd*d + Gw*w ,
%
%  with the Laplace- or Z-transformed plant outputs y, control inputs u, 
%  disturbance input d and noise input w, and with Gu, Gd  and Gw the 
%  corresponding transfer-function matrices.
%  The control input u must correspond to the input group named 'controls'. 
%  If OPTIONS.cdinp = true, then the disturbance input d must correspond 
%  to the input group named 'disturbances'.
%  If OPTIONS.MDFreq is empty or not defined (see below), then DIST is
%  a row vector with N-components, whose j-th entry is computed as follows:
%  DIST(j) = dist(Gu,Gu_j), if OPTIONS.cdinp = false, and
%  DIST(j) = dist([Gu Gd],[Gu_j Gd_j]), if OPTIONS.cdinp = true.
%  The distance function dist(.,.) to be used is defined in accordance with 
%  OPTIONS.distance (see below).
%  FPEAK is a row vector with N components, whose j-th component FPEAK(j) 
%  is the peak frequency (in rad/TimeUnit), where DIST(j) is achieved. 
%  If OPTIONS.MDFreq contains NF frequency values, then DIST is a row
%  vector with N-components, whose j-th entry is the maximum of pointwise 
%  distances between the j-th component models contained in SYSM 
%  and SYS over the NF frequency values. 
%  FPEAK is a row vector with N components, whose j-th component FPEAK(j) 
%  is the peak frequency (in rad/TimeUnit), where DIST(j) is achieved. 
%  MIND is the index of the component model for which the minimum value 
%  of the distances in DIST is achieved. 
%
%  The OPTIONS structure allows to specify various user options, as
%  follows:
%  OPTIONS.tol      - relative tolerance for rank computations
%                     (Default: internally determined value)
%  OPTIONS.distance - specifies the option for the distance function
%                     to be used:
%                     'nugap'    - nu-gap distance (Default) 
%                     'Inf'      - H-infinity norm based distance
%                     '2'        - H2-norm based distance
%                     (Default: 'nugap')
%  OPTIONS.MDFreq   - a set of real frequency values for which the 
%                     pointwise distances have to be computed.
%                     (Default: [])
%  OPTIONS.offset   - the stability boundary offset, to be used to assess  
%                     the finite zeros which belong to the boundary of the 
%                     stability domain as follows: 
%                     in the continuous-time case, these are the finite 
%                        zeros having real parts in the interval 
%                        [-OFFSET, OFFSET]
%                     in the discrete-time case, these are the finite 
%                        zeros having moduli in the interval 
%                        [1-OFFSET,1+OFFSET]. 
%                     (Default: OFFSET = sqrt(eps)).
%  OPTIONS.cdinp    - specifies the option to use both control and 
%                     disturbance input channels to evaluate the 
%                     distances, as follows
%                     true  - use both control and disturbance input channels   
%                     false - use only the control input channels (default)
%
%  See also MDMATCH. 

%  Copyright 2018 A. Varga
%  Author:      A. Varga, 31-10-2018.
%  Revision(s): 
%
%   Method: The evaluation of the nu-gap distances relies on the definition  
%   proposed in [1]. 
%
%   References:
%   [1] G. Vinnicombe. 
%       Uncertainty and feedback: H-infinity loop-shaping and the 
%       nu-gap metric. Imperial College Press, London, 2001. 

narginchk(2,3)
nargoutchk(0,3)

if nargin < 3
    options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

sstype = isa(sysm,'ss');

if ~sstype 
   validateattributes(sysm, {'cell'},{'vector'},'','SYSM')
   N = length(sysm);
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
end

if ~isa(sys,'ss')
   error('SYS must be a LTI state space object')
end

% decode options

% tolerance for rank tests
if isfield(options,'tol')
   tol = options.tol;
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','OPTIONS.tol') 
else
   tol = 0;
end

% zeros selection option
if isfield(options,'distance')
   distance = options.distance;
   validateattributes(distance, {'char'},{'nonempty'},'','OPTIONS.distance') 
else
   distance = 'nugap';
end

switch distance
    case 'nugap'
       jobdist = 0; 
    case {'Inf','inf'}
       jobdist = Inf;
    case '2'
       jobdist = 2;
    otherwise
       error('No such distance selection option')
end


% frequencies for point-wise evaluation of distances
if isfield(options,'MDFreq')
   freq = options.MDFreq;
   if ~isempty(freq)
      validateattributes(freq, {'double'},{'real','vector','>=',0},'','OPTIONS.MDFreq')
   end
else
   freq = [];
end

% offset for the stability region boundary
if isfield(options,'offset')
   offset = options.offset;
   validateattributes(offset, {'double'},{'real','scalar','>',0,'<',1},'','OPTIONS.offset') 
else
   offset = sqrt(eps);
end

% option for extended set of inputs
if isfield(options,'cdinp')
   cdinp = options.cdinp; 
   validateattributes(cdinp, {'logical'},{'binary'},'','OPTIONS.cdinp') 
else
   cdinp = false;
end

if isa(sysm,'ss') 
   [p,~,N] = size(sysm); 
   Ts = sysm.Ts;  
   % decode input information
   if isfield(sysm.InputGroup,'controls')
      % controls
      inpu = sysm.InputGroup.controls;
      m = length(inpu);  
   else
      inpu = []; m = 0;
   end
   tinpu = inpu;
   if cdinp
      if isfield(sysm.InputGroup,'disturbances')
         % disturbances
         inpd = sysm.InputGroup.disturbances;
         m = m + length(inpd);  
         inpu = [inpu inpd];
         tinpd = inpd;
      end
   end
else
   N = length(sysm); 
   p = size(sysm{1},1); 
   Ts = sysm{1}.Ts;  
   for i = 1:N
       % decode input information
       if isfield(sysm{i}.InputGroup,'controls')
          % controls
          tinpu = sysm{i}.InputGroup.controls;
          mu = length(tinpu);  
       else
          tinpu = []; mu = 0;
       end
       if i == 1
          inpu = tinpu; m = mu; 
       else
          if p ~= size(sysm{i},1)
             error('All component models must have the same number of outputs')
          end
          if Ts ~= sysm{i}.Ts;
             error('All component models must have the same sampling time')
          end
          if m ~= mu
              error('All component models must have the same number of control inputs')
          end  
          if ~isequal(inpu,tinpu)
             error('All component models must have the same control inputs')
          end  
       end
   end
   if cdinp
      for i = 1:N
          % decode disturbance input information
          if isfield(sysm{i}.InputGroup,'disturbances')
             % disturbances
             tinpd = sysm{i}.InputGroup.disturbances;
             mdt = length(tinpd);  
          else
             tinpd = []; mdt = 0;
          end
          if i == 1
             inpd = tinpd; md = mdt; m = m+md; 
          else
             if md ~= mdt
              error('All component models must have the same number of disturbance inputs')
             end  
             if ~isequal(inpd,tinpd)
                error('All component models must have the same disturbance inputs')
             end  
          end
      end
      inpu = [inpu inpd];
   end
end

if p ~= size(sys,1)
   error('SYS and all component models in SYSM must have the same number of outputs')
end
if Ts ~= sysm.Ts;
   error('SYS and all component models in SYSM must have the same sampling time')
end
if isfield(sys.InputGroup,'controls')
   % controls
   tu = sys.InputGroup.controls;
else
   tu = [];
end
if ~isequal(tu,tinpu)
   error('SYS and all component models in SYSM must have the same control inputs')
end  
if cdinp
   if isfield(sys.InputGroup,'disturbances')
      td = sys.InputGroup.disturbances;
   else
      td = [];
   end
   if ~isequal(td,tinpd)
      error('SYS and all component models in SYSM must have the same disturbance inputs')
   end  
end

dist = ones(1,N); fpeak = zeros(1,N);
tolinf = max(tol,sqrt(eps)); 
if jobdist 
   for j = 1:N
       if sstype 
          syst = sysm(:,inpu,j)-sys;
       else
          syst = sysm{j}(:,inpu)-sys;
       end
       if isempty(freq)
          % compute the H-inf or H2 norm based distance
          if jobdist == 2
             dist(j) = norm(syst,2);
          else
             [dist(j),fpeak(j)] = norm(syst,inf,tolinf);
          end
       else
          sv = sigma(prescale(syst),freq);
          [dist(j),indfr] = max(sv(1,:),[],2); 
          fpeak(j) = freq(indfr);  
       end
   end
   [~,mind] = min(dist);
   return
end
    
% compute the normalized coprime factorizations of the test system 
% compute the normalized right coprime factorization R1 = [N1;M1]
R1 = grange([sys(:,inpu);eye(m)],struct('inner',true,'tol',tol)); 
% compute the normalized left coprime factorization L1 = [ N1t M1t]
L1 = gcrange([sys(:,inpu) eye(p)],struct('coinner',true,'tol',tol));
for j = 1:N
    if sstype 
       % R2 = [N2;M2]
       R2 = grange([sysm(:,inpu,j);eye(m)],struct('inner',true,'tol',tol));
    else
       R2 = grange([sysm{j}(:,inpu);eye(m)],struct('inner',true,'tol',tol));
    end

    % check conditions on det(R2'*R1)
    syst = gir(R2'*R1,tol);  % syst = R2'*R1 must also work
    [~,infoz] = gzero(syst,tol,offset);
    % check invertibility and presence of zeros on the boundary of 
    % stability domain
    if infoz.nrank ~= order(syst)+m || infoz.nfsbz ;
       continue
    end

    % evaluate winding number 
    [~,infop] = gpole(syst,tol,offset);
    wno = infoz.nfuz - infop.nfuev + infoz.niz - infop.nip;
    % check condition on winding number 
    if wno ~= 0
       % nonzero winding number
       continue
    end
    % compute the underlying system to compute the nu-gap distance 
    % using the definition of Vinnicombe
    syst = L1*[zeros(m,p) -eye(m);eye(p,p+m)]*R2;
    if isempty(freq)
       % compute the nu-gap using the definition of Vinnicombe
       [dist(j),fpeak(j)] = norm(syst,inf,tolinf);
    else
       sv = sigma(prescale(syst),freq);
       [dist(j),indfr] = max(sv(1,:),[],2); 
       fpeak(j) = freq(indfr);  
    end
end

% compute the index of the least distance 
[~,mind] = min(dist);

% end MDDIST2C
