function [mdgain,fpeak,p,relgain]  = mdperf(R,options) 
% MDPERF Model detection distance mapping performance  
%  [MDGAIN,FPEAK] = MDPERF(R,OPTIONS) evaluates for an N x N cell array R 
%  of internal forms of N model detection filters, the M x N array MDGAIN 
%  of gains corresponding to a selected subset of internal forms of 
%  M model detection filters.
%  Each cell R{i,j}, if non-empty, must be given in a standard or 
%  descriptor state-space form, which corresponds to the input-output form
%
%       r_ij = Ru_ij*u + Rd_ij*d_j + Rw_ij*w_j ,
%
%  with the Laplace- or Z-transformed residual outputs r_ij, control inputs
%  u, disturbance inputs d_j and noise inputs w_j, and with Ru_ij, Rd_ij 
%  and Rw_ij the corresponding transfer-function matrices. 
%  The same control inputs u must be defined for all cells of R, 
%  via the input group named 'controls'. If OPTIONS.cdinp = true, then
%  the same disturbance input d is assumed for all component models 
%  (d_j = d) and must correspond to the input group named 'disturbances'. 
%  If OPTIONS.MDFreq is empty or not defined (see below), then 
%  MDGAIN is a M x N non-negative matrix, whose (i,j)-th entry is the 
%  gain of the selected input-output channels of R{k,j}, where k is 
%  the i-th element of OPTIONS.MDSelect. 
%  MDGAIN(i,j) is computed as follows:
%  MDGAIN(i,j) = ||Ru_kj||_inf,         if OPTIONS.cdinp = false, and
%  MDGAIN(i,j) = ||[Ru_kj Rd_kj]||_inf, if OPTIONS.cdinp = true.
%  If OPTIONS.MDFreq contains NF frequency values, then 
%  MDGAIN is an M x N  non-negative matrix, whose (i,j)-th element 
%  MDGAIN(i,j) is computed as the maximum 2-norm gain of the frequency
%  response of the selected input-output channels of R{k,j} over the NF 
%  frequencies contained in FREQ. 
%  FPEAK is an M x N matrix, whose (i,j)-th entry contains the peak 
%  frequency (in rad/TimeUnit), for which MDGAIN(i,j) is achieved. 
%
%  [MDGAIN,FPEAK,P] = MDPERF(R,OPTIONS) evaluates additionally
%  the M x N integer matrix P related to the ordering of gains in 
%  MDGAIN, namely, the i-th row of P is the permutation to be applied to 
%  increasingly reorder the i-th row of MDGAIN. 
%
%  [MDGAIN,FPEAK,P,RELGAIN] = MDPERF(R,OPTIONS) evaluates additionally
%  the M-vector RELGAIN of relative gains. RELGAIN(i) is the ratio 
%  between the second and the l-th smallest gains in the i-th row of 
%  MDGAIN, where l is specified in OPTIONS.MDIndex (see below).
%
%  The OPTIONS structure allows to specify the user options, as follows:
%  OPTIONS.MDSelect - vector of indices of the selected M model detection 
%                     filters for which the gains of the corresponding 
%                     internal forms have to be evaluated 
%                     (Default: 1:N)
%  OPTIONS.MDFreq   - a set of real frequency values for which the maximum  
%                     of 2-norm point-wise gains have to be computed.
%                     (Default: [])
%  OPTIONS.cdinp    - specifies the option to use both control and 
%                     disturbance input channels to evaluate the gains, 
%                     as follows:
%                     true  - use both control and disturbance input channels   
%                     false - use only the control input channels (default)
%  OPTIONS.MDIndex  - specifies the index l of the l-th smallest gains  
%                     to be used to evaluate the relative gains 
%                     to the second smallest gains 
%                     (Default: l = 3)  
%
%  See also MDDIST. 

%  Copyright 2018 A. Varga
%  Author:      A. Varga, 23-10-2018.
%  Revision(s): 
%

narginchk(1,2)
nargoutchk(0,4)

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
if isfield(options,'MDFreq')
   freq = options.MDFreq;
   if ~isempty(freq)
      validateattributes(freq, {'double'},{'real','vector','>=',0},'','OPTIONS.MDFreq')
   end
else
   freq = [];
end

% option to use a specific model to evaluate the relative gains
if isfield(options,'MDIndex')
   MDIndex = options.MDIndex; 
   validateattributes(MDIndex, {'double'},{'integer','>=',1,'<=',N},'','OPTIONS.MDIndex') 
else
   MDIndex = min(3,N);
end

mdgain = zeros(M,N); fpeak = zeros(M,N);
if isempty(freq)
   ki = 0;
   for i = MDSelect
       ki = ki+1;
       for j = 1:N
           if isempty(R{i,j})
              continue
           else
              [mdgain(ki,j),fpeak(ki,j)] = norm(R{i,j}(:,inpu),inf);
           end
       end
   end
else
   ki = 0;
   for i = MDSelect
       ki = ki+1;
       for j = 1:N
           if isempty(R{i,j})
              continue
           else
              sv = sigma(prescale(R{i,j}(:,inpu)),freq); 
              [mdgain(ki,j),indm] = max(sv(1,:),[],2); 
              fpeak(ki,j) = freq(indm);  
           end
       end
   end
end

% compute the gain ordering information
if nargout > 2
    p = zeros(M,N);
    for i = 1:M
       [~,p(i,:)] = sort(mdgain(i,:));
    end
end

% compute the relative gains
if nargout > 3
   relgain = zeros(M,1);
   if N <= 2 || MDIndex == 1, return, end
   for i = 1:M
       relgain(i) = mdgain(i,p(i,2))/mdgain(i,p(i,MDIndex));
   end
end

% end MDPERF
