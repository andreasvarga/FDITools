% Example 5.6 - Solution of an AFDP
% Uses the Control System Toolbox and Descriptor System Tools

% define s as an improper transfer function
s = tf('s');
% define Gu(s), Gw(s), Gf(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     % enter Gu(s)
Gw = [1/(s+2); 0];                   % enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           % enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       % enter dimensions

tol = 1.e-7;   % set tolerance for rank tests

% choose the left nullspace as Q = [I -Gu] and 
% form Rf = Q*[Gf;0] and Rw = Q*[Gw;0]
Q = ss([eye(p) -Gu]); Rf = ss(Gf); Rw = ss(Gw); 

% check solvability using a random frequency
if min(max(abs(evalfr(Rf,rand)))) > 0.01

   % compress Rw to a full row rank matrix
   rw = rank(evalfr(Rw,rand)); nb = size(Q,1);
   if rw < nb
      h = ones(rw,nb);  % usable only for rw = 1
      % use alternatively h = rand(rw,nb);
      Q = h*Q; Rf = h*Rf; Rw = h*Rw;
   end

   % compute the quasi-co-outer-co-inner factorization 
   [Rwi,Rwo]=goifac(Rw,struct('tol',tol)); 

   % compute optimal filter (for standard case)
   Q = gir(Rwo\Q,tol);              % update Q
   Rf = gir(Rwo\Rf,tol); Rw = Rwi;  % update Rf and Rw

   % check for poles on the extended imaginary axis
   poles = gpole([Q Rf]);
   if any(isinf(poles)) || min(abs(real(poles))) < 0.0001

      % compute a stable and proper left coprime factorization
      % of [Q Rf Rw] with desired stability degree -3
      % update if solution is improper or unstable
      opts = struct('sdeg',-3,'smarg',-3); 
      [Q_Rf_Rw,M] = glcf(gir([Q,Rf,Rw],tol),opts);  

      % adjust denominator to unit infinity norm to match example
      Mnorm = -norm(M,inf); 
      Q = minreal(tf(Q_Rf_Rw(:,1:p+mu)/Mnorm),tol)
      Rf = minreal(tf(Q_Rf_Rw(:,p+mu+1:p+mu+mf)/Mnorm),tol)
      Rw = minreal(tf(Q_Rf_Rw(:,p+mu+mf+1:end)/Mnorm),tol)
   end
else
   disp('No solution exists')
end
