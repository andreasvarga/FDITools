function [seli,selord] = ammbasesel(rgain,degs,nout,simple,tol,strongfdi)
% AMMBASESEL Selection of admissible basis vectors to solve the AMMP
% [SELI,SELORD] = AMMBASESEL(RGAIN,DEGS,NOUT,SIMPLE,TOL,STRONGFDI) selects 
% for an NVEC x MR frequency gain matrix RGAIN, several NOUT-tuples of  
% admissible basis vectors to solve the strong FDI problem, if 
% STRONGFDI = TRUE, or the fault detection problem, if STRONGFDI = FALSE. 
% The number of selected vectors must satisfy NOUT <= MIN(NVEC,MR).
% Each row SELI(i,:) contains NOUT indices of basis vectors, such that
% RGAIN(SELI(i,:),:) has full row rank NOUT if STRONGFDI = TRUE, or
% RGAIN(SELI(i,:),:) has all columns nonzero if STRONGFDI = FALSE.
% If the associated NVEC degrees, contained in DEGS, are provided, then
% SELORD(i) is the corresponding tentatively achievable least filter order.
% If SIMPLE = true, a simple basis is assumed, in which case, DEGS(i) is 
% also the order of the i-th basis vector. If SIMPLE = false, a minimum 
% rational basis is assumed. SELORD is empty if DEGS is empty. 
% TOL is a tolerance for rank determinations. 
%

%  Copyright 2018 A. Varga
%  Author:    A. Varga, 17-02-2018.
%  Revisions: A. Varga, 06-08-2018.
%
%  Method: The selection approach in conjunction with the synthesis 
%  Procedure AMMS is described in [1]. 
%
%  References:
%  [1] Varga A.
%      Solving Fault Diagnosis Problems - Linear Synthesis Techniques. 
%      Springer Verlag, 2017.

narginchk(6,6)
nargoutchk(0,2)

validateattributes(rgain, {'double'},{'2d','nonempty'},'','RGAIN',1) 
% nvec  - number of vectors  
% mf    - number of faults  
[nvec,mr] = size(rgain);     

if ~isempty(degs)
   validateattributes(degs(:), {'double'},{'integer','vector','nonnegative','size',[nvec,1]},'','DEGS(:)',2) 
end
% make degs a column vector
if isrow(degs), degs = degs(:); end

validateattributes(strongfdi, {'logical'},{'binary'},'','SIMPLE',6) 
if strongfdi
   validateattributes(nout, {'double'},{'integer','scalar','<=',min(mr,nvec)},'','NOUT',3) 
else
   validateattributes(nout, {'double'},{'integer','scalar','<=',nvec},'','NOUT',3) 
end
validateattributes(simple, {'logical'},{'binary'},'','SIMPLE',4) 
validateattributes(tol, {'double'},{'real','scalar','>=',0},'','TOL',5) 

if nvec == 1
   seli = 1; selord = degs;
   return
end

% find nout combinations of nvec vectors  
seli = nchoosek(1:nvec,nout); ni = size(seli,1);
selord = -ones(ni,1);   % set all orders to -1
ii = false(ni,1); 
for i = 1:ni
    indv = seli(i,:);
    if strongfdi 
        if (tol > 0 && rank(rgain(indv,:),tol) == nout) || ...
           (tol == 0 && rank(rgain(indv,:)) == nout)  
           ii(i) = true;
           if ~isempty(degs)
              % estimate orders 
              selord(i) = sum(degs(indv));
           end
        end
    else
       % evaluate minimum column norm
       beta = norm(rgain(indv,1)); 
       for j = 2:mr
           beta  = min(beta,norm(rgain(indv,j)));
       end  
       if (tol > 0 && beta > tol) || ...
          (tol == 0 && beta > nvec*mr*eps(norm(rgain)))
          ii(i) = true;
          if ~isempty(degs)
             % estimate orders 
             selord(i) = sum(degs(indv));
          end
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
    
% end AMMBASESEL
end

