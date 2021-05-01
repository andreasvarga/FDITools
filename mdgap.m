function [beta,gamma] = mdgap(R,options)
%MDGAP  Noise gaps of model detection filters  
%  GAP = MDGAP(R,OPTIONS) evaluates for an N x N cell array R of internal
%  forms of N model detection filters, the M-dimensional vector GAP 
%  of noise gaps corresponding to a selected subset of internal forms of 
%  M model detection filters. 
%  R{i,j} must be given in a standard or descriptor state-space form, 
%  which corresponds to the input-output form
%
%       r_ij = Ru_ij*u + Rd_ij*d_j + Rw_ij*w_j,
%
%  with the Laplace- or Z-transformed residual outputs r_ij, 
%  control inputs u, disturbance inputs d_j, and noise inputs w_j, and with 
%  Ru_ij, Rd_ij and Rw_ij the corresponding transfer-function matrices.
%  The inputs u, d_j, and w_j of each entry R{i,j} correspond to the  
%  input groups named, respectively, {'controls','disturbances','noise'}.
%  The same control inputs u must be defined for all entries of R, via the  
%  input group named 'controls'. 
%  Note: It is tacitly assumed, that the least infinity norm of the control
%  and disturbance channels of R{i,i}, for i = 1, ...,N, is equal to zero. 
%  GAP(i) is computed as GAP(i) = BETA(i)/GAMMA(i), where for 
%  k representing the i-th element of OPTIONS.MDSelect, 
%  BETA(i) = min_{j~=k}||Ru_kj||_inf,         if OPTIONS.cdinp = false, or
%  BETA(i) = min_{j~=k}||[Ru_kj Rd_kj]||_inf, if OPTIONS.cdinp = true,
%  and GAMMA(i) = ||Rw_kk||_inf.  
%  GAP(i) is set to NaN if the cells R{k,1:N} are empty. 
% 
%  The OPTIONS structure allows to specify the user options, as follows:
%  OPTIONS.MDSelect - vector of indices of the selected M model detection 
%                     filters for which the gaps of the corresponding 
%                     internal forms  have to be evaluated 
%                     (Default: 1:N)
%  OPTIONS.MDFreq   - a set of real frequency values for which the maximum
%                     of 2-norm poit-wise gains have to be computed.
%                     (Default: [])
%  OPTIONS.cdinp    - specifies the option to use both control and 
%                     disturbance input channels to evaluate the gaps, 
%                     as follows:
%                     true  - use both control and disturbance input channels   
%                     false - use only the control input channels (default)
%
%  [BETA,GAMMA] = MDGAP(...) computes explicitly the values of BETA 
%  and GAMMA whose ratio represents the noise gaps. 
%  BETA and GAMMA are M-dimensional vectors and the ratio BETA(i)/GAMA(i) 
%  represents the noise gap of the internal forms k-th filter. 
%  BETA(i) and GAMMA(i) are set to zero if the cells R{k,1:N} are empty. 
%

%  Copyright 2018 A. Varga
%  Author:    A. Varga, 30-10-2018.
%  Revisions: 

narginchk(1,2)
nargoutchk(0,2)

if nargin < 2
    options = struct('freq',[]);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

% option for extended set of inputs
if isfield(options,'cdinp')
   cdinp = options.cdinp; 
   validateattributes(cdinp, {'logical'},{'binary'},'','OPTIONS.cdinp') 
else
   cdinp = false; inpd = [];
end

validateattributes(R, {'cell'},{'2d','square'},'','R')
N = size(R,1);
first = true; 
for i = 1:N
    for j = 1:N
       if isempty(R{i,j})
           continue
       else
           if first
               i1 = i; j1 = j; first = false;
           end
           if ~isa(R{i,j},'ss')
              error('The input R must be cell array of LTI state space objects')
           end
       end
       if i == i1 && j == j1 
          Ts = R{i,j}.Ts; 
          if isfield(R{i,j}.InputGroup,'controls')
             % controls
             inpu = R{i,j}.InputGroup.controls;
          else
             inpu = [];
          end
          if cdinp 
             if isfield(R{i,j}.InputGroup,'disturbances')
                % disturbances
                inpd = R{i,j}.InputGroup.disturbances;
             else
                inpd = [];
             end
          end
       else
          if Ts ~= R{i,j}.Ts
             error('All entries of R must have the same sampling time')
          end
          if isfield(R{i,j}.InputGroup,'controls')
             % controls
             tinpu = R{i,j}.InputGroup.controls;
             if ~isequal(inpu,tinpu)
                error('Different input groups ''control'' defined for two or more entries of R')
             end
          end
          if cdinp 
             if isfield(R{i,j}.InputGroup,'disturbances')
                % disturbances
                tinpd = R{i,j}.InputGroup.disturbances;
                if ~isequal(inpd,tinpd)
                   error('Different input groups ''disturbances'' defined for two or more entries of R')
                end
             end
          end
       end
    end
end    
inpu = [inpu inpd];   

% decode the rest of options

% selected models
if isfield(options,'MDSelect')
   MDSelect = options.MDSelect;
   validateattributes(MDSelect, {'double'},{'vector','integer','>=',1,'<=',N,'increasing'},'','OPTIONS.MDSelect')
else
   MDSelect = 1:N;
end
M = length(MDSelect);

% frequencies for point-wise evaluation of nu-gaps
if isfield(options,'MDfreq')
   freq = options.MDFreq;
   if ~isempty(freq)
      validateattributes(freq, {'double'},{'real','vector','>=',0},'','OPTIONS.MDFreq')
   end
else
   freq = [];
end

if N <= 2
   beta = zeros(M,1); gamma = ones(M,1);
   return
end
   
beta = zeros(M,1); gamma = zeros(M,1);

if isempty(freq)
   ki = 0;
   for i = MDSelect
       ki = ki+1;
       nrm2 = inf;
       for j = 1:N
           if i ~= j
              if isempty(R{i,j})
                 continue
              else
                 nrm2 = min(nrm2,norm(R{i,j}(:,inpu),inf));
              end
           end
       end
       if nrm2 == inf, nrm2 = 0; end
       beta(ki) = nrm2; 
       if isfield(R{i,i}.InputGroup,'noise')
          inpw = R{i,i}.InputGroup.noise; 
          gamma(ki) = norm(R{i,i}(:,inpw),inf);
       end
   end
else
   ki = 0;
   for i = MDSelect
       ki = ki+1;
       nrm2 = inf;
       for j = 1:N
           if i ~= j
              if isempty(R{i,j})
                 continue
              else
                 sv = sigma(prescale(R{i,j}(:,inpu)),freq); 
                 nrm2 = min(nrm2,max(sv(1,:),[],2));
              end
           end
       end
       if nrm2 == inf, nrm2 = 0; end
       beta(ki) = nrm2; 
       if isfield(R{i,i}.InputGroup,'noise')
          inpw = R{i,i}.InputGroup.noise; 
          gamma(ki) = norm(R{i,i}(:,inpw),inf);
       end
   end
end

if nargout <= 1
   beta = beta./gamma;
end
        
% end MDGAP
end
