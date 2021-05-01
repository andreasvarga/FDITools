function [mdgain,fpeak,mind]  = mdmatch(Q,sys,options) 
% MDMATCH Model detection distance matching performance
%  [MDGAIN,FPEAK,MIND] = MDMATCH(Q,SYS,OPTIONS) evaluates for N model 
%  detection filters contained in Q and a given system model SYS, the 
%  N-dimensional vector MDGAIN of gains of the corresponding internal forms.
%  Q is specified as a N x 1 cell array of LTI systems in a standard or 
%  descriptor state-space form, whose i-th component Q{i} corresponds 
%  to the input-output form
%
%       r_i = Q_i*[ y ]
%                 [ u ]
%
%  with the Laplace- or Z-transformed residual outputs r_i, measured 
%  outputs y and control inputs u, and with Q_i the transfer-function 
%  matrix of the i-th filter. The inputs y and  u must correspond to the 
%  input group named 'outputs' and 'controls'. 
%  SYS is LTI system in state-space form, which corresponds to the 
%  input-output form  
%
%      y = Gu*u + Gd*d + Gw*w ,
%
%  with the Laplace- or Z-transformed plant outputs y, control inputs u, 
%  disturbance input d and noise input w, and with Gu, Gd  and 
%  Gw the corresponding transfer-function matrices.
%  The control input u must correspond to the input group named 'controls'. 
%  If OPTIONS.cdinp = true, then the disturbance input d must correspond 
%  to the input group named 'disturbances'.
%  If OPTIONS.MDFreq is empty or not defined (see below), then 
%  MDGAIN is an N-dimensional column vector, whose i-th entry is the 
%  gain of the selected input-output channels of the internal form
%
%        R_i = [ Ru_i Rd_i Rw_i] := Q_i*[ Gu Gd Gw ]
%                                         I  0  0  ]
%  MDGAIN(i) is computed as follows:
%  MDGAIN(i) = ||Ru_i||_inf,         if OPTIONS.cdinp = false, and
%  MDGAIN(i) = ||[Ru_i Rd_i]||_inf,  if OPTIONS.cdinp = true.
%  If OPTIONS.MDFreq contains NF frequency values, then 
%  MDGAIN is an N-dimensional column vector, whose i-th entry MDGAIN(i) 
%  is computed as the maximum 2-norm gain of the frequency
%  response of the selected input-output channels of R_i over the NF 
%  frequencies contained in FREQ. 
%  FPEAK is an an N-dimensional column vector, whose i-th entry contains 
%  the peak frequency (in rad/TimeUnit), for which MDGAIN(i,j) is achieved. 
%  MIND is the index of the component of MDGAIN for which the minimum value
%  of the gains is achieved. For a properly designed filter Q, this is
%  also the index of the best matching model to the current model SYS. 
%
%  The OPTIONS structure allows to specify the user options, as follows:
%  OPTIONS.MDFreq   - a set of real frequency values for which the 
%                     maximum of 2-norm gains have to be computed.
%                     (Default: [])
%  OPTIONS.cdinp    - specifies the option to use both control and 
%                     disturbance input channels to evaluate the gains, 
%                     as follows:
%                     true  - use both control and disturbance input channels   
%                     false - use only the control input channels (default)
%
%  See also MDDIST2C. 

%  Copyright 2018 A. Varga
%  Author:      A. Varga, 31-10-2018.
%  Revision(s): 
%

narginchk(2,3)
nargoutchk(0,3)

if nargin < 3
    options = struct('freq',[]);
else
   validateattributes(options,{'struct'},{'nonempty'},'','OPTIONS')
end

% option for extended set of inputs
if isfield(options,'cdinp')
   cdinp = options.cdinp; 
   validateattributes(cdinp, {'logical'},{'binary'},'','OPTIONS.cdinp') 
else
   cdinp = false; sinpd = [];
end


validateattributes(Q, {'cell'},{'vector'},'','Q')
N = length(Q);
first = true; 
for i = 1:N
    if isempty(Q{i})
       continue
    else
       if first
          i1 = i;  first = false;
       end
       if ~isa(Q{i},'ss')
          error('The input Q must be cell array of LTI state space objects')
       end
    end
    if i == i1  
       Ts = Q{i}.Ts; 
       if isfield(Q{i}.InputGroup,'controls')
          % controls
          inpu = Q{i}.InputGroup.controls;
       else
          inpu = [];
       end
       if isfield(Q{i}.InputGroup,'outputs')
          % outputs
          inpy = Q{i}.InputGroup.outputs;
       else
          inpy = [];
       end
    else
       if Ts ~= Q{i}.Ts
          error('All entries of Q must have the same sampling time')
       end
       if isfield(Q{i}.InputGroup,'controls')
          % controls
          tinpu = Q{i}.InputGroup.controls;
          if ~isequal(inpu,tinpu)
             error('Different input groups ''controls'' defined for two or more entries of Q')
          end
       end
       if isfield(Q{i}.InputGroup,'outputs')
          % controls
          tinpy = Q{i}.InputGroup.outputs;
          if ~isequal(inpy,tinpy)
             error('Different input groups ''outputs'' defined for two or more entries of Q')
          end
       end
    end
end    

if ~isa(sys,'ss')
   error('SYS must be a LTI state space object')
end

if length(inpy) ~= size(sys,1) 
   error('The number of outputs of SYS and of inputs in the input groups ''outputs'' of Q must be the same')
end
if Ts ~= sys.Ts;
   error('SYS and the component filters of Q must have the same sampling time')
end
if isfield(sys.InputGroup,'controls')
   % controls
   sinpu = sys.InputGroup.controls;
else
   sinpu = [];
end
if length(inpu) ~= length(sinpu)
   error('The number of control imputs of SYS and of inputs in the input groups ''controls'' of Q must be the same')
end  
if cdinp
   if isfield(sys.InputGroup,'disturbances')
      sinpd = sys.InputGroup.disturbances;
   end
end
mu = length(sinpu);
md = length(sinpd);

% decode options

% frequencies for point-wise evaluation of nu-gaps
if isfield(options,'MDFreq')
   freq = options.MDFreq;
   if ~isempty(freq)
      validateattributes(freq, {'double'},{'real','vector','>=',0},'','OPTIONS.MDFreq')
   end
else
   freq = [];
end

mdgain = zeros(N,1); fpeak = zeros(N,1);
inpQ = [inpy inpu]; syse = [sys(:,[sinpu sinpd]); eye(mu,mu+md)];
if isempty(freq)
   for i = 1:N
       if isempty(Q{i})
          continue
       else
          [mdgain(i),fpeak(i)] = norm(Q{i}(:,inpQ)*syse,inf);
       end
   end
else
   for i = 1:N
       if isempty(Q{i})
          continue
       else
          sv = sigma(prescale(Q{i}(:,inpQ)*syse),freq); 
          [mdgain(i),indm] = max(sv(1,:),[],2); 
          fpeak(i) = freq(indm);  
       end
   end
end

% compute the index of best matching model
[~,mind] = min(mdgain);

% end MDMATCH
