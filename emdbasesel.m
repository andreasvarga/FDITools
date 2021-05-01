function [seli,selord] = emdbasesel(S,degs,rdim,nout,simple)
% EMDBASESEL Selection of admissible basis vectors for solving the EMDP
% [SELI,SELORD] = EMDBASESEL(S,DEGS,RDIM,NOUT,SIMPLE) selects several 
% NOUT-tuples of admissible basis vectors, using the NVEC x LFREQ x N 
% binary structure matrix S, such that the exact model detection problem 
% (EMDP) is solvable by using model detection filters with RDIM outputs.
% S is a three-dimensional matrix, such that, for each k, S(:,:,k) has
% nonzero columns. LFREQ > 1 is employed, when several structure 
% vectors are provided at different frequency values 
% (see the output of MDSSPEC).   
% Each row SELI(i,:) contains NOUT indices of basis vectors, whose
% linear combination is admissible, i.e. , S(SELI(i,:),:,k) has all columns
% nonzero for all k. 
% If the associated NVEC degrees contained in DEGS, are provided, then
% SELORD(i) is the corresponding tentatively achievable least filter order.
% If SIMPLE = true, a simple basis is assumed, in which case, DEGS(i) is 
% also the order of the i-th basis vector. If SIMPLE = false, a minimum 
% rational basis is assumed.
% SELORD is empty if DEGS is empty. 
%
% 
%  Copyright 2017-2018 A. Varga
%  Author:     A. Varga, 26-08-2017.
%  Revisions:  
%
%  Method: The selection approach, used in conjunction with the synthesis 
%  Procedure EMD, is described in [1]. 
%
%  References:
%  [1] Varga A.
%      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%      Springer Verlag, 2017.

narginchk(5,5)
nargoutchk(0,2)

validateattributes(S, {'logical'},{'binary','3d','nonempty'},'','S',1) 
[nvec,lfreq,N] = size(S);   % number of vectors  

if ~isempty(degs)
   validateattributes(degs(:), {'double'},{'integer','vector','nonnegative','size',[nvec,1]},'','DEGS(:)',2) 
end
% make degs a column vector
if isrow(degs), degs = degs(:); end

validateattributes(rdim, {'double'},{'integer','scalar','>=',1,'<=',nvec},'','RDIM',3)   
validateattributes(nout, {'double'},{'integer','scalar','>=',rdim,'<=',nvec},'','NOUT',4)   
validateattributes(simple, {'logical'},{'binary'},'','SIMPLE',5) 

% quick return if possible
if nvec == 1
   seli = 1; selord = degs;
   return
end

% find rdim combinations of nout vectors which solve the EFDP 
seli = nchoosek(1:nvec,nout); ni = size(seli,1);
selord = -ones(ni,1);   % set all orders to -1
nqmax = sum(degs);
ii = true(ni,1); 
for i = 1:ni
    indv = seli(i,:);
    % check admissibility
    stest = true;
    for j = 1:N
        for k = 1:lfreq
            stest = stest & any(S(indv,k,j));
        end
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
    
% end EMDBASESEL
end

