function [beta,gamma] = fdifscond(R,freq,S)
%FDIFSCOND  Fault sensitivity condition of FDI filters 
%  FSCOND = FDIFSCOND(R) computes FSCOND, the fault sensitivity condition 
%  of the internal form of a fault detection filter R or the vector FSCOND 
%  of fault sensitivity conditions of the internal forms of a collection 
%  of fault detection and isolation filters.
%  If R is a stable LTI filter, R has the input-output form
%
%       r = Rf*f + Rv*v ,                            (1)
%
%  where r, f and v are the Laplace- or Z-transformed residual vector, 
%  fault input vector, and the additional (non-relevant) inputs v, 
%  respectively, and where Rf(lambda) and Rv(lambda) are the 
%  transfer function matrices corresponding to the respective inputs.
%  The input f of R corresponds to the input group named 'faults'. 
%  FSCOND is computed as FSCON = BETA/GAMMA, where BETA is the 
%  H-(infinity minus)-index of Rf (i.e., the least H-infinity norm
%  of columns) and GAMMA is the maximum H-infinity norm of columns of Rf. 
%  FSCOND = 0 if one or more fault inputs are undetectable. 
%
%  If R is an N x 1 array of stable LTI systems, R{i}, i = 1, ..., N, where 
%  the i-th system R{i} has an input-output form
%
%       y_i = Rf_i*f + Rv_i*v ,                            (2)
%  
%  then FSCOND is an N-dimensional vector, with FSCOND(i) representing the 
%  fault sensitivity condition of the i-th system R{i}. FSCOND(i) is 
%  computed as FSCOND(i) = BETA(i)/GAMMA(i), where BETA(i) is the 
%  H-(infinity minus)-index of Rf_i and GAMMA(i) is the maximum of the 
%  H-infinity norm of the columns of Rf_i. 
%
%  FSCOND = FDIFSCOND(R,FREQ) computes, for the stable LTI system R  
%  in (1), the fault sensitivity condition as FSCOND = BETA/GAMMA, where 
%  BETA is the H-(infinity minus)-index of Rf evaluated over the  
%  frequencies contained in the real vector FREQ, and GAMMA is the GAMMA is
%  the maximum of the 2-norm of columns of Rf evaluated over the  
%  frequencies contained in the real vector FREQ. 
% 
%  If R is an N x 1 array of LTI systems, R{i}, i = 1, ..., N, with 
%  the i-th system R{i} having the input-output form in (2), 
%  then FSCOND is an N-dimensional vector, with FSCOND(i) representing the 
%  fault sensitivity condition of the i-th system R{i}. FSCOND(i) is 
%  computed as FSCOND(i) = BETA(i)/GAMMA(i), where BETA(i) is the 
%  H-(infinity minus)-index of Rf_i evaluated over the frequencies 
%  contained in the real vector FREQ and GAMMA(i) is the maximum of 
%  the 2-norms of the columns of Rf_i, evaluated over the frequencies 
%  contained in the real vector FREQ. FSCOND(i) is set to NaN if R{i} 
%  is empty. 
%
%  FSCOND = FDIFSCOND(R,[],S) computes, for the stable LTI system R  
%  in (1) with q outputs and mf fault inputs, and a q x mf 
%  logical structure matrix S, the q-dimensional vector FSCOND, whose i-th 
%  element FSCOND(i) is the fault sensitivity condition computed 
%  (for the i-th rows) as FSCOND(i) = BETA(i)/GAMMA(i), with BETA(i),  
%  the H-(infinity minus)-index of Rf1_i and GAMMA(i), the maximum of  
%  the H-infinity norm of the columns of Rf1_i, where Rf1_i is formed from 
%  the elements of the i-th row of Rf for which S(i,j) = true.
%
%  If R is an N x 1 array of stable LTI systems, with R{i} as in (2)
%  and S is an N x mf logical structure matrix, then FSCOND is an 
%  N-dimensional vector, whose i-th element FSCOND(i) is the fault 
%  sensitivity condition of the i-th system R{i}, computed as 
%  FSCOND(i) = BETA(i)/GAMMA(i), with BETA(i), the H-(infinity minus)-index
%  of Rf1_i and GAMMA(i) the maximum of the H-infinity norm of the columns 
%  of Rf1_i, where Rf1_i is formed from the columns of Rf_i for which 
%  S(i,j) = true. FSCOND(i) is set to NaN if R{i} is empty. 
%  
%  FSCOND = FDIFSCOND(R,FREQ,S) computes, for the stable LTI system R
%  in (1) with p outputs and mf fault inputs, and a p x mf 
%  logical structure matrix S, the p-dimensional vector FSCOND, whose i-th 
%  element FSCOND(i) is the fault sensitivity condition (for the i-th rows)
%  computed as FSCOND(i) = BETA(i)/GAMMA(i), with BETA(i), 
%  the H-(infinity minus)-index of Rf1_i evaluated over the frequencies
%  contained in the real vector  FREQ,  and GAMMA(i) the maximum of the 
%  2-norms of the columns of Rf1_i, evaluated over the frequencies 
%  contained in the real vector FREQ, where Rf1_i is formed from the 
%  elements of the i-th row of Rf for which S(i,j) = true.
%
%  If R is an N x 1 array of stable LTI systems, with R{i} as in (2)
%  and S is an N x mf logical structure matrix, then FSCOND is an 
%  N-dimensional vector, whose i-th element FSCOND(i) is the fault 
%  sensitivity condition of the i-th system R{i}, computed as 
%  FSCOND(i) = BETA(i)/GAMMA(i), with BETA(i), the H-(infinity minus)-index
%  of Rf1_i evaluated over the frequencies contained in the real vector 
%  FREQ,  and GAMMA(i), the maximum of the 2-norms of the columns of Rf1_i,
%  evaluated over the frequencies contained in the real vector FREQ, where 
%  Rf1_i is formed from the columns of Rf_i for which S(i,j) = true. 
%  FSCOND(i) is set to NaN if R{i} is empty. 
%  
%  [BETA,GAMMA] = FDIFSCOND(...) computes explicitly the values of BETA 
%  and GAMMA whose ratio represents the fault sensitivity condition. 
%  If R is an N x 1 array of LTI systems, then BETA and GAMMA
%  are N-dimensional vectors and the ratio BETA(i)/GAMA(i) represents
%  the fault sensitivity condition of the i-th component systems R{i}. 
%  BETA(i) and GAMMA(i) are set to zero if R{i} is empty. 

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
       if ~isa(R{i},'ss')
          error('R must be cell array of state-space LTI system objects')
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
   end
elseif isa(R,'ss')
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
else
   error('R must be either a cell array of LTI systems or a LTI system object')
end

if mf  == 0
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
      gamma = hinfmax(R(:,inpf),freq);
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
          beta(i) = hinfminus(R(i,inpfi),freq);
          gamma(i) = hinfmax(R(i,inpfi),freq);
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
          beta(i) = hinfminus(R{i},freq);
          gamma(i) = hinfmax(R{i},freq);
       else
          sysc = R{i};
          inpfi = inpf(S(i,:));
          beta(i) = hinfminus(sysc(:,inpfi),freq);
          gamma(i) = hinfmax(sysc(:,inpfi),freq);
       end
   end
   if nargout <= 1
      beta = beta./gamma;
   end
end

% end FDIFSCOND
end
