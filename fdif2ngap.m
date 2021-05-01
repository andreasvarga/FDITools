function [beta,gamma] = fdif2ngap(R,freq,S)
%FDIF2NGAP  Fault-to-noise gap of FDI filters  
%  GAP = FDIF2NGAP(R) computes GAP, the fault-to-noise gap of the internal 
%  form of a fault detection filter R or the vector GAP of 
%  fault-to-noise gaps of the internal forms of a collection of 
%  fault detection and isolation filters.
%  If R is a stable LTI filter, R has the input-output form
%
%       r = Rf*f + Rw*w + Rv*v ,                            (1)
%
%  with the Laplace- or Z-transformed residual outputs r, fault inputs f, 
%  noise inputs w, and auxiliary inputs v, and with Rf, Rw and Rv the  
%  corresponding transfer function matrices. The inputs f and w of R 
%  correspond to the input groups named 'faults' and 'noise', respectively. 
%  GAP is computed as GAP = BETA/GAMMA, where BETA is the 
%  H-(infinity minus)-index of Rf (i.e., the least H-infinity norm
%  of columns) and GAMMA is the H-infinity norm of Rw. GAP is infinite if
%  there are no noise inputs and GAP = 0 if there are no fault inputs. 
%
%  If R is an N x 1 array of stable LTI systems, R{i}, i = 1, ..., N, where 
%  the i-th system R{i} has an input-output form
%
%       y_i = Rf_i*f + Rw_i*w ,                    (2)
%
%  then GAP is an N-dimensional vector, with GAP(i) representing the 
%  fault-to-noise gap of the i-th system R{i}. GAP(i) is computed as
%  GAP(i) = BETA(i)/GAMMA(i), where BETA(i) is the H-(infinity minus)-index
%  of Rf_i and GAMMA(i) is the H-infinity norm of Rw_i. If R{i} is
%  empty, then GAP(i) is set to NaN.
%
%  GAP = FDIF2NGAP(R,[],S) computes, for the stable LTI system R in 
%  (1) with q outputs, mf fault inputs and mw noise inputs, and a q x mf 
%  logical structure matrix S, the q-dimensional vector GAP, whose i-th 
%  element GAP(i) is the fault-to-noise gap computed (for the i-th rows) as 
%  GAP(i) = BETA(i)/GAMMA(i), with BETA(i), the H-(infinity minus)-index
%  of Rf1(i,:) and GAMMA(i), the H-infinity norm of [Rf2(i,:) Rw(i,:)], 
%  where Rf1 and Rf2 are formed from the columns of Rf for which 
%  S(i,j) = true and, respectively, S(i,j) = false.
%
%  If R is an N x 1 array of stable LTI systems, with R{i} as in (2)
%  and S is an N x mf logical structure matrix, then GAP is an N-dimensional 
%  vector, whose i-th element GAP(i) is the fault-to-noise gap of 
%  the i-th system R{i}, computed as GAP(i) = BETA(i)/GAMMA(i), with 
%  BETA(i), the H-(infinity minus)-index of Rf1_i and GAMMA(i) the 
%  H-infinity norm of [ Rf2_i Rw_i], where Rf1_i and Rf2_i are formed from 
%  the columns of Rf_i for which S(i,j) = true and, respectively, 
%  S(i,j) = false. If R{i} is empty, then GAP(i) is set to NaN.
%  
%  GAP = FDIF2NGAP(R,FREQ) computes, for the stable LTI system R in 
%  (1), the fault-to-noise gap as GAP = BETA/GAMMA, where BETA is the 
%  H-(infinity minus)-index of Rf evaluated over the frequencies contained  
%  in the real vector FREQ, and GAMMA is the H-infinity norm of Rw.
% 
%  If R is an N x 1 array of LTI systems, R{i}, i = 1, ..., N, with 
%  the i-th system R{i} having the input-output form in (2), 
%  then GAP is an N-dimensional vector, with GAP(i) representing the 
%  fault-to-noise gap of the i-th system R{i}. GAP(i) is computed as
%  GAP(i) = BETA(i)/GAMMA(i), where BETA(i) is the H-(infinity minus)-index
%  of Rf_i evaluated over the frequencies contained in the real vector FREQ
%  and GAMMA(i) is the H-infinity norm of Rw_i. If R{i} is empty, then 
%  GAP(i) is set to NaN.
%
%  GAP = FDIF2NGAP(R,FREQ,S) computes, for the stable LTI system R in 
%  (1) with q outputs, mf fault inputs and mw noise inputs, and a q x mf 
%  logical structure matrix S, the q-dimensional vector GAP, whose i-th 
%  element GAP(i) is the fault-to-noise gap computed (for the i-th rows) as 
%  GAP(i) = BETA(i)/GAMMA(i), with BETA(i), the H-(infinity minus)-index
%  of Rf1(i,:) evaluated over the frequencies contained in the real vector 
%  FREQ,  and GAMMA(i), the H-infinity norm of [Rf2(i,:) Rw(i,:)], 
%  where Rf1 and Rf2 are formed from the columns of Rf for which 
%  S(i,j) = true and, respectively, S(i,j) = false.
%
%  If R is an N x 1 array of stable LTI systems, with R{i} as in (2)
%  and S is an N x mf logical structure matrix, then GAP is an N-dimensional 
%  vector, whose i-th element GAP(i) is the fault-to-noise gap of 
%  the i-th system R{i}, computed as GAP(i) = BETA(i)/GAMMA(i), with 
%  BETA(i), the H-(infinity minus)-index of Rf1_i evaluated over the 
%  frequencies contained in the real vector FREQ,  and GAMMA(i) the 
%  H-infinity norm of [ Rf2_i Rw_i], where Rf1_i and Rf2_i are formed from 
%  the columns of Rf_i for which S(i,j) = true and, respectively, 
%  S(i,j) = false. If R{i} is empty, then GAP(i) is set to NaN.
%  
%  [BETA,GAMMA] = FDIF2NGAP(...) computes explicitly the values of BETA 
%  and GAMMA whose ratio represents the fault-to-noise gap. 
%  If R is an N x 1 array of LTI systems, then BETA and GAMMA
%  are N-dimensional vectors and the ratio BETA(i)/GAMA(i) represents
%  the fault-to-noise gap of the i-th component systems R{i}. 
%  If R{i} is empty, then BETA(i) and GAMMA(i) are set to zero.

%  Copyright 2018 A. Varga
%  Author:    A. Varga, 06-07-2018.
%  Revisions: A. Varga, 01-08-2018, 16-10-2018.

narginchk(1,3)
nargoutchk(0,2)

% check input system form
if isa(R,'cell')
   N = length(R);
   validateattributes(R, {'cell'},{'vector'},'','R')
   % select the indices of non-empty systems 
   syssel = zeros(1,N);
   for i = 1:N
       if ~isempty(R{i}) 
          syssel(i) = i;
       end 
   end
   syssel = syssel(syssel > 0); 
   if isempty(syssel)
      error('R must be a non-empty cell array of LTI system objects')
   end
   for i = syssel
       if ~isa(R{i},'lti')
          error('R must be cell array of LTI system objects')
       end 
   end
   j = syssel(1); 
   Ts = R{j}.Ts;
   if isfield(R{1}.InputGroup,'faults')
      % faults
      inpf = R{j}.InputGroup.faults;
      mf = length(inpf);
   else
      mf = size(R{j},2); 
      inpf = 1:mf; 
   end    
   if isfield(R{j}.InputGroup,'noise')
      % faults
      inpw = R{j}.InputGroup.noise;
      mw = length(inpw);
   else
      mw = 0; 
      inpw = []; 
   end    
   for i = syssel(2:end)
       if Ts ~= R{i}.Ts
          error('All component models must have the same sampling time')
       end
       if isfield(R{1}.InputGroup,'faults')
          % fault inputs
          inpfi = R{i}.InputGroup.faults;
          mfi = length(inpfi);
       else
          mfi = size(R{i},2); 
          inpfi = 1:mfi; 
       end    
       if mf ~= mfi || ~isequal(inpf,inpfi)
          error('All component models must have the same number of fault inputs')
       end
       if isfield(R{1}.InputGroup,'noise')
          % noise inputs
          inpwi = R{i}.InputGroup.noise;
          mwi = length(inpwi);
       else
          mwi = 0; 
          inpwi = []; 
       end    
       if mw ~= mwi || ~isequal(inpw,inpwi)
          error('All component models must have the same number of noise inputs')
       end
   end
elseif isa(R,'lti')
   [~,~,N] = size(R);
   if N ~= 1
       error('No multiple models supported yet')
   end
   if isfield(R.InputGroup,'faults')
      % faults
      inpf = R.InputGroup.faults;
      mf = length(inpf);
   else
      mf = 0;
      inpf = []; 
   end    
   if isfield(R.InputGroup,'noise')
      % faults
      inpw = R.InputGroup.noise;
      mw = length(inpw);
   else
      mw = 0;
      inpw = []; 
   end    
else
   error('R must be either a cell array or a LTI system object')
end

m = mf+mw;
if m == 0
   if N == 1 
      beta = []; gamma = [];
   else
     beta = zeros(N,0); gamma = zeros(N,0);
   end
   return
end

if nargin < 2
    freq = [];
else
    if ~isempty(freq)
        validateattributes(freq, {'double'},{'real','vector','>=',0},'','FREQ',2)
    end
end

if nargin < 3
   S = [];
else
   if ~isempty(S)
       validateattributes(S, {'logical'},{'binary'},'','S',3)
       if mf ~= size(S,2)
          error('Structure matrix incompatible with the number of faults')
       end
       nb = size(S,1);
       if N == 1
          if size(R,1) ~= nb 
             error('Structure matrix incompatible with the number of outputs')
          end
       else
          if N ~= nb 
             error('Structure matrix incompatible with the number of systems')
          end
       end
   end
end

if isa(R,'lti')
   if isempty(S) 
      beta = hinfminus(R(:,inpf),freq);
      gamma = norm(R(:,inpw),inf);
      if isinf(gamma)
          error('The system R has a pole on the boundary of the stability region')
      end
      if nargout <= 1
         beta = beta/gamma;
      end
   else
      beta = zeros(nb,1); gamma = zeros(nb,1);
      for i = 1:nb
          inpfi = inpf(S(i,:));
          inpwi = inpf(~S(i,:));
          beta(i) = hinfminus(R(i,inpfi),freq);
          gamma(i) = norm(R(i,[inpwi inpw]),inf);
      end
      if any(isinf(gamma))
          error('The system R has a pole on the boundary of the stability region')
      end
      if nargout <= 1
         beta = beta./gamma;
      end
   end
   return
elseif isa(R,'cell')
   beta = zeros(N,1); gamma = zeros(N,1);    
   for i = syssel   
       if isempty(S)
          [beta(i),gamma(i)] = fdif2ngap(R{i},freq);
       else
          sysc = R{i};
          sysc.InputGroup.faults = inpf(S(i,:));
          sysc.InputGroup.noise = [inpf(~S(i,:)) inpw];         
          [beta(i),gamma(i)] = fdif2ngap(sysc,freq);
       end
   end
   if nargout <= 1
      beta = beta./gamma;
   end
end

% end FDIF2NGAP
end
