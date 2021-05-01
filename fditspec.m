function smat = fditspec(R,tol,FDTol,freq,blkopt) 
%FDITSPEC Computation of the weak or strong structure matrix 
%  SMAT = FDITSPEC(R,TOL,FDTOL) determines for a LTI system R,
%  representing the internal form of a fault detection and isolation (FDI) 
%  filter, the binary weak structure matrix SMAT corresponding to the  
%  nonzero elements of the transfer function matrix Rf in the 
%  input-output form of R
%
%      r = Rf*f + Rv*v                                (1)
%
%  where r, f and v are the Laplace- or Z-transformed q-dimensional 
%  residual vector, mf-dimensional fault input vector, and the 
%  additional (non-relevant) inputs v, respectively, and where Rf(lambda)
%  and Rv(lambda) are the transfer-function matrices corresponding 
%  to the respective inputs.
%  The inputs f of R must correspond to the input groups named 'faults'.
%  SMAT is determined as a q x mf logical matrix, whose (i,j)-th element 
%  is SMAT(i,j) = true, if the (i,j)-th element of Rf(lambda) is 
%  nonzero, and otherwise, SMAT(i,j) = false. R is assumed to be in a 
%  descriptor system representation with the realization 
%  (Ar-lambda*Er,Bfr,Cr,Dfr) of Rf(lambda). Otherwise, an automatic 
%  conversion to a descriptor realization is internally performed. 
%  TOL is a relative tolerance used for controllability tests (a default 
%  value is internally computed if TOL = 0). FDTOL is an absolute  
%  threshold for the magnitudes of nonzero elements in the system matrices 
%  Bfr, Cr and Dfr. If FDTOL = 0, the default value 
%  FDTOL = 0.0001*max([1,norm(Bfr,1),norm(Cr,inf),norm(Dfr,1)]) is used. 
%  If the input group 'faults' is not defined for the system R, then SMAT  
%  is the structure matrix corresponding to transfer function matrix of R.
%
%  SMAT = FDITSPEC(R,TOL,FDTOL,FREQ) determines for a LTI system R 
%  representing the internal form of a FDI filter and nf real frequency 
%  values contained in FREQ, the binary strong structure matrix SMAT 
%  corresponding to the nonzero elements of the transfer function matrix
%  Rf(lambda) in (1) for all frequency values contained in the vector FREQ. 
%  SMAT is determined as a q x mf x nf logical matrix, whose k-th page 
%  corresponds to the nonzero elements of the transfer-function matrix 
%  Rf(lambda) evaluated in the k-th frequency value FREQ(k). 
%  Accordingly, SMAT(i,j,k) = true, if the (i,j)-th element of 
%  Rf(lambda) has no zero in the complex frequency lambda(k) corresponding 
%  to the real frequency value FREQ(k) (see Note below). 
%  Otherwise, SMAT(i,j,k) = false. 
%  R is assumed to be in a descriptor system representation with the
%  realization (Ar-lambda*Er,Bfr,Cr,Dfr) of Rf(lambda).  Otherwise, an 
%  automatic conversion to a descriptor realization is internally performed. 
%  TOL is a relative tolerance used for observability and controllability 
%  tests (a default value is internally computed if TOL = 0). 
%  FDTOL is an absolute threshold for the magnitudes of nonzero elements 
%  in the system matrices Bfr, Cr and Dfr, and also for the rank tests on 
%  the system matrix [Ar-lambda*Er Bfr; Cr Dfr]. 
%  If FDTOL = 0, the following default values are used: 
%  FDTOL = 0.0001*max([norm(Bfr,1),norm(Cr,inf),norm(Dfr,1)]) 
%  is used for the magnitude of nonzero elements in Bfr, Cr, and Dfr. 
%  For the rank tests of the system matrix, the value 
%  FDTOL = 0.0001*max(norm(1,[Ar Bfr;Cr Dfr],1),norm(Er,1)) is used. 
%  If the input group 'faults' is not defined for the system R, then SMAT  
%  is the structure matrix corresponding to transfer function matrix of R.
%
%  SMAT = FDITSPEC(R,TOL,FDTOL,[],'block') determines for a LTI system 
%  R, representing the internal form of a FDI filter, the binary weak 
%  structure (row) vector SMAT of the transfer function matrix Rf(lambda)
%  in (1) using a block-structure  based evaluation.  
%  SMAT is determined as a 1 x mf row vector, whose j-th element 
%  SMAT(1,j) = true if the j-th column of Rf(lambda) is  nonzero, and 
%  SMAT(1,j) = false if the j-th column of Rf(lambda) is zero. 
%  If the input group 'faults' is not defined for the system R, then SMAT  
%  is the structure matrix corresponding to transfer function matrix of R.
%
%  SMAT = FDITSPEC(R,TOL,FDTOL,FREQ,'block') determines, for a LTI system 
%  R, representing the internal form of a FDI filter, and for nf real 
%  frequency values contained in FREQ, the binary (strong) structure (row)  
%  vector SMAT using a block-wise analysis of the transfer function matrix 
%  Rf(lambda) in (1) for all frequency values contained in the vector FREQ. 
%  SMAT is determined as a 1 x mf x nf logical matrix, whose k-th page 
%  corresponds to the nonzero columns of the transfer-function 
%  matrix Rf(lambda) evaluated in the k-th frequency value FREQ(k). 
%  Accordingly, SMAT(1,j,k) = true if the j-th column of 
%  Rf(lambda) has no zero in the complex frequency corresponding to the real
%  frequency value FREQ(k) (see Note below). Otherwise, SMAT(1,j,k) = false.   
%  If the input group 'faults' is not defined for the system R, then SMAT  
%  is the structure matrix corresponding to transfer function matrix of R.
%
%  SMAT = FDITSPEC(R,TOL,FDTOL) determines for a cell array of LTI 
%  systems R{i}, i = 1, ..., N, the binary (weak) structure matrix SMAT, 
%  whose i-th row is the structure vector which results by using a 
%  block-wise analysis of zero/nonzero columns of G_i(lambda), the 
%  transfer function matrix of R{i}. 
%  Accordingly, SMAT(i,j) = true if the j-th column of G_i(lambda) is   
%  nonzero and SMAT(i,j) = false if the j-th column of G_i(lambda) is zero.   
%  The LTI systems R{i}, i = 1, ..., N, must have the same number of 
%  inputs and the same sampling time. 
%  If the same input group 'faults' is defined for all component systems 
%  R{i}, i = 1, ..., N, then SMAT is the structure matrix corresponding 
%  to R{i}(:,'faults'), for i = 1, ..., N. 
%  All entries of the i-th row of SMAT are set to false if R{i} is empty. 
%
%  SMAT = FDITSPEC(R,TOL,FDTOL,FREQ) determines for a cell array of LTI 
%  systems R{i}, i = 1, ..., N, representing the internal form of 
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
%  Accordingly, SMAT(i,j,k) = true if the j-th column of Rf_i(lambda) has 
%  no zero in the complex frequency corresponding to the real frequency 
%  value FREQ(k) (see Note below). Otherwise, SMAT(i,j,k) = false.  
%  If the same input group 'faults' is defined for all component systems 
%  R{i}, i = 1, ..., N, then SMAT is the structure matrix corresponding 
%  to R{i}(:,'faults'), for i = 1, ..., N. 
%  If the input group 'faults' is not defined for the system R{i}, then 
%  the i-th row of SMAT is the structure vector corresponding to 
%  transfer function matrix of R{i}.
%  If R{i} is empty, then all entries of the i-th page SMAT(i,:,:) are set  
%  to false. 
%
%  Note: The complex frequency corresponding to a real frquency value 
%  FREQ(k) is lambda(k) := sqrt(-1)*FREQ(k), for a continuous-time system 
%  R, and lambda(k) := exp(sqrt(-1)*FREQ(k)*abs(R.Ts)) for a discrete-time
%  system R. 
%
%  See also FDISSPEC.

%  Copyright 2016-2018 A. Varga
%  Author:    A. Varga, 02-12-2016.
%  Revisions: A. Varga, 09-02-2017, 21-08-2017, 06-07-2018, 24-10-2018,
%                       15-12-2018, 12-06-2019.
%
%  Method: For the definition of the structure matrix, see [1]. For the
%  determination of the weak structure matrix, controllable realizations
%  are determined for each column of G(lambda) and the nonzero elements 
%  in each column are identified  (see Corollary 7.1 of [1]).
%  For the determination of the strong structure matrix, minimal
%  realizations are determined for each element of G(lambda) and checks on
%  the absence of zeros are performed by checking the full rank of the
%  system matrix for all frequencies in FREQ (see Corollary 7.2 in [1]).
%  For the block-structure based determination of the weak structure 
%  matrix, the non-zero columns of G(lambda) are identified from
%  their controllable realizations (see Corollary 7.1 of [1]), while for
%  the block-structure based determination of the strong structure 
%  matrix checks on the absence of zeros of the columns of G(lambda) 
%  are performed by checking the full rank of the system matrix 
%  of their controllable realizations for all frequencies in FREQ 
%  (see Corollary 7.2 in [1]).
%
%  References:
%  [1] Varga A.
%      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%      Springer Verlag, 2017; sec.3.4.

narginchk(1,5)
nargoutchk(0,1)

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
   validateattributes(tol, {'double'},{'real','scalar','>=',0},'','TOL',2) 
else
   tol = 0;
end

if nargin > 2
   validateattributes(FDTol, {'double'},{'real','scalar','>=',0},'','FDTOL',3) 
else
   FDTol = 0;
end

if nargin > 3
   if ~isempty(freq) 
      validateattributes(freq, {'double'},{'real','vector','>=',0},'','FREQ',4)
   end
else
   freq = [];
end

if nargin > 4
   if ~isempty(blkopt) 
      validateattributes(blkopt, {'char'},{},'','BLKOPT',5)
   end
   block = strcmp(blkopt,'block');
   if ~block
      error('Improper option: use instead ''block'' to specify block-wise analysis')
   end
else
   block = false;
end

if nargin < 2 || tol < 0
    tol = 0;
end

if nargin < 3 
    FDTol = 0;
end

if nargin < 4
   freq = [];
end


if isa(R,'lti')
   % convert to state-space if necessary
   if ~isa(R,'ss')
       R = gir(ss(R),tol);   
   end
   % set strong analysis option
   Sstrong = ~isempty(freq);
   if Sstrong
      lfreq = length(freq);
      w = 1i*freq;                   % w = j*freq
      Ts = abs(R.Ts);
      if Ts > 0, w = exp(Ts*w); end  % w = exp(j*Ts*freq)
   end   
   for isys = 1:N 
       % perform scaling to equilibrate norms of system matrices
       warning('off','all')
       % [a,b,c,d,e] = dssdata(ssbal(sys(:,inpf,isys,1)));
       [a,b,c,d,e] = dssdata(R(:,inpf,isys,1));
       warning('on','all')
       % elliminate possible non-dynamic modes 
       [a,e,b,c,d] = sl_gminr(2,a,e,b,c,d,tol,0);  
       if FDTol > 0
          FDSTol = FDTol;
       else
          FDTol = 0.0001*max([1, norm(b,1),norm(c,inf),norm(d,1)]);
          if Sstrong 
             FDSTol = 0.0001*max([1, norm([a b; c d],1), norm(e,1)]);
          end
       end
       
       swmat =  abs(d) > FDTol; 
       if Sstrong
          [p,mf]=size(d); 
          if block
             smat = false(1,mf,lfreq);
          else
             smat = false(p,mf,lfreq);
          end
       else
          smat = swmat; 
          if isempty(a) ||  all(smat(:)) || (block && all(max(smat,[],1)))
             % finish for a constant system or all entries/columns nonzero
             if block
                smat = max(smat,[],1);
             end
             return
          end
       end
   
       % employ structural analysis to compute weak/strong structure matrix  
       iaout = find(max(abs(c),[],2)> FDTol)'; % select nonzero rows in C
       jainp = find(max(abs(b),[],1)> FDTol);  % select nonzero columns in B
       for j = jainp 
           % elliminate uncontrollable eigenvalues for the i-th column of B
           [a1,e1,b1,c1,d1] = sl_gminr(1,a,e,b(:,j),c,d(:,j),tol,1,0);
           swmat(iaout,j) = abs(d1(iaout,:)) > FDTol | ...
                            max(abs(c1(iaout,:)),[],2) > FDTol;
           if Sstrong
              if block
                 for k = 1:lfreq
                     if any(swmat(:,j))
                        % check if freq(k) is a zero of the (i,j)-th element 
                        if isinf(freq(k))
                           % check for infinite zero 
                           s = svd([e1,b1;c1 d1]);
                        else
                           % check for a finite zero 
                           s = svd([a1-w(k)*e1,b1;c1 d1]);
                        end
                        smat(1,j,k) = s(end) > FDSTol;
                     else
                        break 
                     end
                 end
                  
              else
                 isel = find(swmat(:,j)); 
                 for i = isel 
                     % elliminate unobservable eigenvalues for the i-th row of C
                     [a2,e2,b2,c2,d2] = sl_gminr(1,a1,e1,b1,c1(i,:),d1(i,:),tol,2,0);
                     for k = 1:lfreq
                         if swmat(i,j)
                            % check if freq(k) is a zero of the (i,j)-th element 
                            if isinf(freq(k))
                               % check for infinite zero 
                               s = svd([e2,b2;c2 d2]);
                            else
                               % check for a finite zero 
                               s = svd([a2-w(k)*e2,b2;c2 d2]);
                            end
                            smat(i,j,k) = s(end) > FDSTol;
                         else
                            break 
                         end
                     end
                 end
              end
           end
       end
       if ~Sstrong
          if block
             smat = max(swmat,[],1);
          else
             smat = swmat;
          end
       end
   end
elseif isa(R,'cell')
   lfreq = max(1,length(freq)); 
   smat = false(N,mf,lfreq);
   for i = syssel
       smat(i,:,:) = fditspec(R{i}(:,inpf),tol,FDTol,freq,'block');
   end
end
   
% end FDITSPEC
end
