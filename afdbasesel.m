function [seli,selord] = afdbasesel(S,rwgain,degs,rdim,nout,simple,tol)
% AFDBASESEL Selection of admissible basis vectors for solving the AFDP
% [SELI,SELORD] = AFDBASESEL(S,RWGAIN,DEGS,RDIM,NOUT,SIMPLE,TOL) selects 
% several NOUT-tuples of admissible basis vectors, using the binary 
% structure matrix S and the frequency gain matrix RWGAIN,
% such that the approximate fault detection problem (AFDP) is 
% solvable by using fault detection filters with RDIM outputs.
% If rank(RWGAIN) > 0, then RDIM must not exceed rank(RWGAIN).
% S is a NVEC x MF x N (3-dimensional) matrix, such that, for each k, 
% S(:,:,k) has all columns nonzero. MF is the number of faults and  N > 1 
% is employed, when several structure matrices are provided at different  
% frequency values (see the output of FDISSPEC).  
% RWGAIN is a NVEC x MW matrix, where MW is the number of noise inputs.   
% Each row SELI(i,:) contains NOUT indices of basis vectors, whose
% linear combination is admissible, i.e. , S(SELI(i,:),:,k) has all columns
% nonzero for all k and, if RWGAIN is nonzero, then  RWGAIN(SELI(i,:),:) 
% has full row rank.   
% If the associated NVEC degrees contained in DEGS are provided, then
% SELORD(i) is the corresponding tentatively achievable least filter order.
% If SIMPLE = true, a simple basis is assumed, in which case, DEGS(i) is 
% also the order of the i-th basis vector. If SIMPLE = false, a minimum 
% rational basis is assumed.
% SELORD is empty if DEGS is empty. 
% TOL is a non-negative tolerance used for rank determinations. 
%
% See also FDISSPEC.

% Copyright 2018 A. Varga
% Author:     A. Varga, 30-04-2018.
% Revisions:  
%
% Method: The selection approach, used in conjunction with the synthesis 
% Procedure AFD, is described in [1]. 
%
% References:
% [1] Varga A.
%     Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%     Springer Verlag, 2017.

narginchk(7,7)
nargoutchk(0,2)

%validateattributes(S, {'logical'},{'binary','3d','nonempty'},'','S',1) 
validateattributes(S, {'logical'},{'binary','3d'},'','S',1) 
[nvec,~,n] = size(S);   % numbers of vectors, faults, structure matrices  
% nvec  - number of vectors  
% mf    - number of faults 
% n     - number of structure matrices

if isempty(rwgain)
   mw = 0;
else
   validateattributes(rwgain, {'double'},{'2d','nonempty'},'','RWGAIN',2) 
   [nvec1,mw] = size(rwgain); % mw is the number of noise inputs  
   if nvec ~= nvec1 
       error('S and RWGAIN must have the same number of rows')
   end
end

if ~isempty(degs)
   validateattributes(degs(:), {'double'},{'integer','vector','nonnegative','size',[nvec,1]},'','DEGS(:)',3) 
end
% make degs a column vector
if isrow(degs), degs = degs(:); end

validateattributes(rdim, {'double'},{'integer','scalar','>=',1,'<=',nvec},'','RDIM',4)   
validateattributes(nout, {'double'},{'integer','scalar','>=',rdim,'<=',nvec},'','NOUT',5)   
validateattributes(simple, {'logical'},{'binary'},'','SIMPLE',6) 
validateattributes(tol, {'double'},{'real','scalar','>=',0},'','TOL',7) 

% determine rank of RWGAIN
if mw 
   if tol > 0,
       rw = rank(rwgain,tol);
   else
       rw = rank(rwgain);
   end
   if rw && rdim > rw
       error('RDIM must not exceed the rank of RWGAIN')
   end       
else
    rw = 0;
end

if nvec == 1 && rw == 0
   seli = 1; selord = degs;
   return
end


% find rdim combinations of nout vectors which solve the AFDP 
seli = nchoosek(1:nvec,nout); ni = size(seli,1);
selord = -ones(ni,1);   % set all orders to -1
nqmax = sum(degs);
ii = true(ni,1); 
for i = 1:ni
    indv = seli(i,:);
    % check admissibility
    stest = true;
    for k = 1:n
        stest = stest & all(max(S(indv,:,k),[],1));
    end
    if stest
       if ~isempty(degs)
          % estimate orders 
          if simple || rdim == nout
             % degree = the sums of degrees of selected vectors
             selord(i) = sum(degs(indv));
          else
             % degree = RDIM times the maximum degree of selected vectors
             selord(i) = min(nqmax,rdim*max(degs(indv)));
          end
       end
    else
       ii(i) = false;
    end
    % check full rank conditions
    if rw && ii(i)
       if (tol > 0 && rank(rwgain(indv,:),tol) < rdim) || ...
          (tol == 0 && rank(rwgain(indv,:)) < rdim)     
          ii(i) = false;
       end
    end
end

seli = seli(ii,:);

if isempty(degs)
   selord = [];
else
   selord = selord(ii);
   % sort row combinations to ensure increasing tentative orders  
   [selord,ii] = sort(selord);
   seli = seli(ii,:);
end
    
% end AFDBASESEL
end
