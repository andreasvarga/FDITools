function [smat,gains] = fdisspec(R,FDGainTol,freq,blkopt) 
%FDISSPEC Computation of the strong structure matrix 
%  [SMAT,GAINS] = FDISSPEC(R,FDGAINTOL,FREQ) determines for a LTI 
%  system R representing the internal form of a fault detection and 
%  isolation (FDI) filter and nf real frequency values contained in FREQ,
%  the binary strong structure matrix SMAT corresponding to the fault
%  inputs in the input-output form of R
%
%      r = Rf*f + Rv*v                                (1)
%
%  where r, f and v are the Laplace- or Z-transformed q-dimensional 
%  residual vector, mf-dimensional fault input vector, and the 
%  additional (non-relevant) inputs v, respecetively, and where Rf(lambda)
%  and Rv(lambda) are the transfer-function matrices corresponding 
%  to the respective inputs.
%  The inputs f of R must correspond to the input groups named 'faults'.
%  SMAT is determined as a q x mf x nf logical matrix, whose k-th page 
%  corresponds to the nonzero elements of the transfer-function matrix 
%  Rf(lambda) evaluated in the k-th frequency value FREQ(k). 
%  Accordingly, SMAT(i,j,k) = true, if the frequency-response gain of the 
%  (i,j)-th element of Rf(lambda(k)), is greater than the positive threshold
%  FDGAINTOL, where lambda(k) is the complex frequency corresponding 
%  to the real frequency value FREQ(k). Otherwise, SMAT(i,j,k) = false.
%  For a continuous-time system R, lambda(k) = sqrt(-1)*FREQ(k), and for 
%  a discrete-time system R, lambda(k) = exp(sqrt(-1)*FREQ(k)*abs(R.Ts)). 
%  The frequency values in FREQ must be such that none of the corresponding 
%  complex frequency values is a pole of Rf. 
%  GAINS is a q x mf non-negative matrix, whose (i,j)-th element contains 
%  the minimum value of the frequency-response gains of the (i,j)-th 
%  element of Rf(lambda) evaluated over all frequencies in FREQ 
%  (also known as the H-minus index of the (i,j)-th element).
%  If the input group 'faults' is not defined for the system R, then SMAT  
%  is the structure matrix corresponding to transfer function matrix of R.
%
%  [SMAT,GAINS] = FDISSPEC(R,FDGAINTOL,FREQ,'block') determines for a LTI 
%  system R representing the internal form of a FDI filter and nf real 
%  frequency values contained in FREQ, the binary strong structure matrix 
%  SMAT corresponding a block-wise check of the presence of zeros of the 
%  columns of Rf(lambda) in (1) in the frequency values contained in FREQ.  
%  SMAT is determined as a 1 x mf x nf logical matrix, whose k-th page 
%  corresponds to the nonzero columns of the transfer-function 
%  matrix Rf(lambda) evaluated in the k-th frequency value FREQ(k). 
%  Accordingly, SMAT(1,j,k) = true, if the norm of the frequency-response 
%  of the j-th column of Rf(lambda(k)), is greater than the positive 
%  threshold FDGAINTOL, where lambda(k) is the complex frequency 
%  corresponding  to the real frequency value FREQ(k) (see above). 
%  Otherwise, SMAT(1,j,k) = false.
%  GAINS is a non-negative mf-dimensional row vector, whose j-th element 
%  contains the minimum norm of the frequency-responses of the j-th column 
%  of Rf(lambda) evaluated over all frequencies in FREQ (also known as the 
%  H-minus index of the j-th column).
%  If the input group 'faults' is not defined for the system R, then SMAT  
%  is the structure matrix corresponding to transfer function matrix of R.
%
%  [SMAT,GAINS] = FDISSPEC(R,FDGAINTOL,FREQ) determines for a cell array  
%  of LTI systems R{i}, i = 1, ..., N, representing the internal form of 
%  a collection of N FDI filters and nf real frequency values 
%  contained in FREQ, the binary structure matrix SMAT, whose i-th row is 
%  the structure vector corresponding to the fault inputs in the 
%  input-output form of R{i}
%
%      r_i = Rf_i*f + Rv_i*v                        (2)        
%
%  where Rf_i(lambda) and Rv_i(lambda) are the transfer-function matrices 
%  corresponding to the respective inputs. The inputs f of R{i} must 
%  correspond to the input groups named 'faults', and all LTI systems 
%  R{i}, i = 1, ..., N, must have the same number of fault inputs and 
%  the same sampling time. 
%  The i-th row of SMAT is the structure vector which results by using a 
%  block-wise check of the presence of zeros of the columns of Rf_i(lambda),
%  in the frequency values contained in FREQ. 
%  Accordingly, SMAT(i,j,k) = true, if the norm of the frequency-response 
%  of the j-th column of Rf_i(lambda(k)), is greater than the positive 
%  threshold FDGAINTOL, where lambda(k) is the complex frequency 
%  corresponding  to the real frequency value FREQ(k) (see above). 
%  Otherwise, SMAT(i,j,k) = false.
%  If the input group 'faults' is not defined for the system R{i}, then 
%  the i-th row of SMAT is the structure vector corresponding to 
%  transfer function matrix of R{i}.
%  If R{i} is empty, then all entries of the i-th page SMAT(i,:,:) are set  
%  to false. 
%  GAINS is a N x mf matrix, whose (i,j)-th entry contains the minimum of 
%  he norms of the frequency-responses of the j-th column of Rf_i evaluated
%  over all complex frequencies corresponding to FREQ (see above). 
%  If R{i} is empty, then all entries of the i-th row of GAINS are set to
%  zero.
%
%  [SMAT,GAINS] = FDISSPEC(R) uses the default values 
%  FDGAINTOL = 0.01 and FREQ = 0 to determine SMAT corresponding to the 
%  magnitudes of the DC-gains of Rf(lambda) in (1) or Rf_i in (2). 
%
%  See also FDITSPEC.

%  Copyright 2016-2018 A. Varga
%  Author:    A. Varga, 02-12-2016.
%  Revisions: A. Varga, 26-08-2017, 10-05-2018, 06-07-2018, 19-10-2018,
%                       15-12-2018, 12-06-2019.
%
%  Method: For the definitions of the structure matrix and H-minus index,
%  see [1].  
%
%  References:
%  [1] Varga A.
%      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%      Springer Verlag, 2017; sec.3.4.

narginchk(1,4)
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
      smat = false(N,0);
      return
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
      m = size(R{j},2); 
      inpf = 1:m; 
      mf = m;
   end    
   for i = syssel(2:end)
       if Ts ~= R{i}.Ts
          error('All component models must have the same sampling time')
       end
       if isfield(R{1}.InputGroup,'faults')
          % faults
          inpfi = R{i}.InputGroup.faults;
          mfi = length(inpfi);
       else
          mi = size(R{i},2); 
          inpfi = 1:mi;
          mfi = mi;
       end    
       if mf ~= mfi || ~isequal(inpf,inpfi)
          error('All component models must have the same number of fault inputs')
       end
   end
elseif isa(R,'ss') || isa(R,'tf') || isa(R,'zpk')
   [~,m,N] = size(R);
   if N ~= 1
       error('No multiple models supported yet')
   end
   if isfield(R.InputGroup,'faults')
      % faults
      inpf = R.InputGroup.faults;
      mf = length(inpf);
   else
      inpf = 1:m; 
      mf = m;
   end    
else
   error('R must be either a cell array or a LTI system object')
end

if nargin > 1
   validateattributes(FDGainTol, {'double'},{'real','scalar','>=',0},'','FDGAINTOL',2) 
   if FDGainTol == 0
      FDGainTol = 0.01;
   end    
else
   FDGainTol = 0.01;
end

if nargin > 2
   validateattributes(freq, {'double'},{'real','vector','>=',0},'','FREQ',3)
else
   freq = 0;
end

if nargin > 3
   if ~isempty(blkopt) 
      validateattributes(blkopt, {'char'},{},'','BLKOPT',4)
   end
   block = strcmp(blkopt,'block');
   if ~block
      error('Improper option: use instead ''block'' to specify block-wise analysis')
   end
else
   block = false;
end

if isa(R,'lti')
   p = size(R(:,inpf,1),1); 
   lfreq = length(freq);
   if mf == 0 
      if block
         smat = false(1,0);
         gains = zeros(1,0);
      else
         smat = false(p,0);
         gains = zeros(p,0);
      end
      return
   end
   
   Ts = abs(R.Ts);
   discr = (Ts > 0);
   if discr
      w = exp(Ts*1i*freq);   % w = exp(j*Ts*freq)
   else
      w = 1i*freq;           % w = j*freq
   end  
   for isys = 1:N 

       % compute strong structure matrix of R by evaluating the element-wise 
       % H- index over the given frequency values

       gs = evalfr(R(:,inpf,isys,1),w(1));
       if any(isinf(gs(:)))
          error('FDISSPEC:pole',['The frequency ',num2str(w(1)),' is a system pole'])
       end
       if block
          smat = false(1,mf,lfreq); 
          gains = zeros(1,mf);
          for j = 1:mf 
              gains(1,j) = norm(gs(:,j));
          end
          smat(1,:,1) = (gains > FDGainTol);
          for i = 2:lfreq
              gs = evalfr(R(:,inpf,isys,1),w(i));
              if any(isinf(gs(:)))
                 error('FDISSPEC:pole',['The frequency ',num2str(w(i)),' is a system pole'])
              end
              for j = 1:mf 
                  nrmgsj = norm(gs(:,j)); 
                  smat(1,j,i) = (nrmgsj > FDGainTol);
                  gains(1,j) = min(gains(1,j),nrmgsj);
              end
          end
       else
          smat = false(p,mf,lfreq); 
          gains = abs(gs);
          smat(:,:,1) = (gains > FDGainTol);
          for i = 2:lfreq
              gs = abs(evalfr(R(:,inpf,isys,1),w(i)));
              if any(isinf(gs(:)))
                 error('FDISSPEC:pole',['The frequency ',num2str(w(i)),' is a system pole'])
              end
              smat(:,:,i) = (gs > FDGainTol);
              ind = find(gs < gains);
              gains(ind) = gs(ind);
          end
       end
   end
elseif isa(R,'cell')
   lfreq = max(1,length(freq)); 
   smat = false(N,mf,lfreq);
   gains = zeros(N,mf);
   for i = syssel
       [smat(i,:,:),gains(i,:)] = fdisspec(R{i}(:,inpf),FDGainTol,freq,'block');
   end
end

% end FDISSPEC
end
