function [dist,fpeak,perm,reldist]  = mddist(sysm,options) 
% MDDIST Computation of distances between component models 
%  [DIST,FPEAK] = MDDIST(SYSM,OPTIONS) evaluates for a multiple 
%  model SYSM with N component models the distances DIST between the 
%  N component models and a subset of M <= N component models. 
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
%  If OPTIONS.MDFreq is empty or not defined (see below), then DIST is a
%  M x N non-negative matrix, whose (i,j)-th entry is the distance 
%  between the selected input-output channels of the k-th and j-th 
%  component models contained in SYSM, where k is the i-th element of  
%  OPTIONS.MDSelect. 
%  DIST(i,j) is computed as follows:
%  DIST(i,j) = dist(Gu_k,Gu_j), if OPTIONS.cdinp = false, and
%  DIST(i,j) = dist([Gu_k Gd_k],[Gu_j Gd_j]), if OPTIONS.cdinp = true.
%  The distance function dist(.,.) to be used is defined in accordance with 
%  OPTIONS.distance (see below).
%  FPEAK is a M x N matrix, whose (i,j)-th entry contains the peak 
%  frequency (in rad/TimeUnit), where NUGAP(i,j) is achieved. 
%  If OPTIONS.MDFreq contains NF frequency values, then 
%  DIST is a M x N  non-negative matrix, whose (i,j)-th element
%  DIST(i,j) is the maximum of pointwise distances between the 
%  selected channels of the  k-th and j-th component models contained in 
%  SYSM over the NF frequency values. 
%  FPEAK is a M x N matrix, whose (i,j)-th entry contains the peak 
%  frequency (in rad/TimeUnit), for which the pointwise distance 
%  achieves its peak value. 
%
%  [DIST,FPEAK,PERM] = MDDIST(SYSM,OPTIONS) evaluates additionally
%  the M x N integer matrix PERM related to the ordering of distances in 
%  DIST, namely, the i-th row of PERM is the permutation to be applied to 
%  increasingly reorder the i-th row of DIST.
%
%  [DIST,FPEAK,PERM,RELDIST] = MDDIST(SYSM,OPTIONS) evaluates additionally
%  the M-vector RELDIST of relative distances, such that,  RELDIST(i) is 
%  the ratio of the second and the l-th smallest distances in the 
%  i-th row of DIST, where l is specified in OPTIONS.MDIndex (see below).
%
%  The OPTIONS structure allows to specify various user options, as
%  follows:
%  OPTIONS.MDSelect - vector of indices of the selected M models to which 
%                     the distances have to be evaluated 
%                     (Default: 1:N)
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
%  OPTIONS.MDIndex  - specifies the index l of the l-th smallest  
%                     distances to be used to evaluate the relative 
%                     distances to the second smallest distances 
%                     (Default: l = 3)  
%  OPTIONS.cdinp    - specifies the option to use both control and 
%                     disturbance input channels to evaluate the 
%                     distances, as follows
%                     true  - use both control and disturbance input channels   
%                     false - use only the control input channels (default)
%
%  See also MDPERF. 

%  Copyright 2018 A. Varga
%  Author:      A. Varga, 31-10-2018.
%  Revision(s): 
%
%   Method: The evaluation of the nu-gap distances relies on the definition  
%   proposed in [1]. For efficiency purposes, the intervening normalized
%   factorizations of the components systems are performed only once and
%   all existing symmetries are exploited.
%
%   References:
%   [1] G. Vinnicombe. 
%       Uncertainty and feedback: H-infinity loop-shaping and the 
%       nu-gap metric. Imperial College Press, London, 2001. 

narginchk(1,2)
nargoutchk(0,4)

if nargin < 2
   options = struct('tol',0);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

sstype = isa(sysm,'ss');

if sstype 
   [~,~,N] = size(sysm);
elseif isa(sysm,'cell')
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
else
   error('SYSM must be either a cell array or a LTI state space object')
end

% decode options

% selected models
if isfield(options,'MDSelect')
   MDSelect = options.MDSelect;
   validateattributes(MDSelect, {'double'},{'vector','integer','>=',1,'<=',N,'increasing'},'','OPTIONS.MDSelect')
else
   MDSelect = 1:N;
end
M = length(MDSelect);

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

% option to use a specific model to evaluate the relative distances
if isfield(options,'MDIndex')
   MDIndex = options.MDIndex; 
   validateattributes(MDIndex, {'double'},{'integer','>=',1,'<=',N},'','OPTIONS.MDIndex') 
else
   MDIndex = min(3,N);
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
   % decode input information
   if isfield(sysm.InputGroup,'controls')
      % controls
      inpu = sysm.InputGroup.controls;
      m = length(inpu);  
   else
      inpu = []; m = 0;
   end
   if cdinp
      if isfield(sysm.InputGroup,'disturbances')
         % disturbances
         inpd = sysm.InputGroup.disturbances;
         m = m + length(inpd);  
         inpu = [inpu inpd];
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

dist = zeros(M,N); 
fpeak = zeros(M,N);
tolinf = max(tol,sqrt(eps)); 
indN = 1:N;
if jobdist 
   ki = 0;
   for i = MDSelect
       ki = ki+1;
       jind = setdiff(indN,MDSelect(1:ki));
       ji = intersect(MDSelect(ki+1:end),jind);
       for j = jind
           if j ~= i
              if sstype 
                 syst = sysm(:,inpu,j)-sysm(:,inpu,i);
              else
                 syst = sysm{j}(:,inpu)-sysm{i}(:,inpu);
              end
              if isempty(freq)
                 % compute the H-inf or H2 norm based distance
                 if jobdist == 2
                    dist(ki,j) = norm(syst,2);
                 else
                    [dist(ki,j),fpeak(ki,j)] = norm(syst,inf,tolinf);
                 end
              else
                 sv = sigma(prescale(syst),freq);
                 [dist(ki,j),indfr] = max(sv(1,:),[],2); 
                 fpeak(ki,j) = freq(indfr);  
              end
           end
       end
       if ~isempty(ji)
          dist(ki+1:end,i) = dist(ki,ji)';
          fpeak(ki+1:end,i) = fpeak(ki,ji)';
       end
   end
   % compute the distance ordering information 
   if nargout > 2
      perm = zeros(M,N);
      for i = 1:M
          [~,perm(i,:)] = sort(dist(i,:));
      end
   end

   % compute the relative distances
   if nargout > 3
      reldist = zeros(M,1);
      if N <= 2 || MDIndex == 1, return, end
      for i = 1:M
          reldist(i) = dist(i,perm(i,2))/dist(i,perm(i,MDIndex));
      end
   end
   return
end

% compute the normalized right coprime factorizations only once
R2 = cell(1,N);
for j = 1:N
    if sstype 
       % R2 = [N2;M2]
       R2{j} = grange([sysm(:,inpu,j),;eye(m)],struct('inner',true,'tol',tol));
    else
       R2{j} = grange([sysm{j}(:,inpu),;eye(m)],struct('inner',true,'tol',tol));
    end
end

ki = 0;
for i = MDSelect
    ki = ki+1;
    % compute normalized left and right factorizations of the i-th system 
    if sstype 
       % compute the normalized right coprime factorization R1 = [N1;M1]
       R1 = grange([sysm(:,inpu,i);eye(m)],struct('inner',true,'tol',tol)); 
       % compute the normalized left coprime factorization L1 = [ N1t M1t]
       L1 = gcrange([sysm(:,inpu,i) eye(p)],struct('coinner',true,'tol',tol));
    else
       % compute the normalized right coprime factorization R1 = [N1;M1]
       R1 = grange([sysm{i}(:,inpu);eye(m)],struct('inner',true,'tol',tol)); 
       % compute the normalized left coprime factorization L1 = [ N1t M1t]
       L1 = gcrange([sysm{i}(:,inpu) eye(p)],struct('coinner',true,'tol',tol));
    end
    jind = setdiff(indN,MDSelect(1:ki));
    ji = intersect(MDSelect(ki+1:end),jind);
    for j = jind
        if j ~= i
%            if sstype 
%               % R2 = [N2;M2]
%               R2 = grange([sysm(:,inpu,j),;eye(m)],struct('inner',true,'tol',tol));
%            else
%               R2 = grange([sysm{j}(:,inpu),;eye(m)],struct('inner',true,'tol',tol));
%            end

           % check conditions on det(R2'*R1)
           syst = gir(R2{j}'*R1,tol);  % syst = R2'*R1 must also work
           [~,infoz] = gzero(syst,tol,offset);
           % check invertibility and presence of zeros on the boundary of 
           % stability domain
           if infoz.nrank ~= order(syst)+m || infoz.nfsbz ;
              dist(ki,j) = 1; fpeak(ki,j) = 0;
              continue
           end

           % evaluate winding number 
           [~,infop] = gpole(syst,tol,offset);
           wno = infoz.nfuz - infop.nfuev + infoz.niz - infop.nip;
           % check condition on winding number 
           if wno ~= 0
              % nonzero winding number
              dist(ki,j) = 1; fpeak(ki,j) = 0;
              continue
           end
           % compute the underlying system to compute the nu-gap distance 
           % using the definition of Vinnicombe
           syst = L1*[zeros(m,p) -eye(m);eye(p,p+m)]*R2{j};
           if isempty(freq)
              % compute the nu-gap using the definition of Vinnicombe
              [dist(ki,j),fpeak(ki,j)] = norm(syst,inf,tolinf);
           else
              sv = sigma(prescale(syst),freq);
              [dist(ki,j),indm] = max(sv(1,:),[],2); 
              fpeak(ki,j) = freq(indm);  
           end
        end
    end
    if ~isempty(ji)
       dist(ki+1:end,i) = dist(ki,ji)';
       fpeak(ki+1:end,i) = fpeak(ki,ji)';
    end
end

% compute the distance ordering information 
if nargout > 2
    perm = zeros(M,N);
    for i = 1:M
       [~,perm(i,:)] = sort(dist(i,:));
    end
end

% compute the relative distances
if nargout > 3
   reldist = zeros(M,1);
   if N <= 2 || MDIndex == 1, return, end
   for i = 1:M
       reldist(i) = dist(i,perm(i,2))/dist(i,perm(i,MDIndex));
   end
end
       

% end MDDIST
