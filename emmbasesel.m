function [seli,selord] = emmbasesel(rgain,degs,nout,simple,tol)
% EMMBASESEL Selection of admissible basis vectors to solve the strong EFDIP
% [SELI,SELORD] = EMMBASESEL(RGAIN,DEGS,NOUT,SIMPLE,TOL) selects for an 
% NVEC x MF frequency gain matrix RGAIN with full column rank, several 
% NOUT-tuples of admissible basis vectors, such that the strong fault 
% detection and isolation problem (strong EFDIP) is solvable by using 
% fault detection filters with MF outputs. 
% The number of selected vectors must satisfy NOUT >= MF.
% Each row SELI(i,:) contains NOUT indices of basis vectors, such that
% RGAIN(SELI(i,:),:) is full column rank. 
% If the associated NVEC degrees contained in DEGS, are provided, then
% SELORD(i) is the corresponding tentatively achievable least filter order.
% If SIMPLE = true, a simple basis is assumed, in which case, DEGS(i) is 
% also the order of the i-th basis vector. If SIMPLE = false, a minimum 
% rational basis is assumed. SELORD is empty if DEGS is empty. 
% TOL is a tolerance for rank determinations. 
%

%  Copyright 2017-2018 A. Varga
%  Author:    A. Varga, 05-04-2017.
%  Revisions: A. Varga, 25-08-2017.
%
%  Method: The selection approach in conjunction with the synthesis 
%  Procedure EMMS is described in [1]. 
%
%  References:
%  [1] Varga A.
%      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%      Springer Verlag, 2017.

narginchk(5,5)
nargoutchk(0,2)

validateattributes(rgain, {'double'},{'2d','nonempty'},'','RGAIN',1) 
% nvec  - number of vectors  
% mf    - number of faults  
[nvec,mf] = size(rgain);     

if ~isempty(degs)
   validateattributes(degs(:), {'double'},{'integer','vector','nonnegative','size',[nvec,1]},'','DEGS(:)',2) 
end
% make degs a column vector
if isrow(degs), degs = degs(:); end

validateattributes(nout, {'double'},{'integer','scalar','>=',mf,'<=',nvec},'','NOUT',3)   
validateattributes(simple, {'logical'},{'binary'},'','SIMPLE',4) 
validateattributes(tol, {'double'},{'real','scalar','>=',0},'','TOL',5) 

if nvec == 1
   seli = 1; selord = degs;
   return
end

% find rdim combinations of nout vectors which solve the strong EFDIP 
seli = nchoosek(1:nvec,nout); ni = size(seli,1);
selord = -ones(ni,1);   % set all orders to -1
nqmax = sum(degs);
ii = true(ni,1); 
for i = 1:ni
    indv = seli(i,:);
    if (tol > 0 && rank(rgain(indv,:),tol) == mf) || ...
       (tol == 0 && rank(rgain(indv,:)) == mf)     
       if ~isempty(degs)
          % estimate orders 
          if simple || mf == nout
             % degree = the sums of degrees of selected vectors
             selord(i) = sum(degs(indv));
          else
             % degree = MF times the maximum degree of selected vectors
             selord(i) = min(nqmax,mf*max(degs(indv)));
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
    
% end EMMBASESEL
end

