% Example 5.11 - Solution of an AFDIP

% Uses the Control System Toolbox and Descriptor System Tools (DSTOOLS)

% define s as an improper transfer function
s = tf('s');
% define Gu(s), Gw(s), Gf(s)
Gu = [(s+1)/(s+2); (s+2)/(s+3)];     % enter Gu(s)
Gw = [1/(s+2); 0];                   % enter Gw(s)
Gf = [(s+1)/(s+2) 0; 0 1];           % enter Gf(s)
p = 2; mu = 1; mw = 1; mf = 2;       % enter dimensions
S = eye(mf);                         % enter structure matrix

% Procedure AFDI

% Step 1): choose the left nullspace as Q1 = [I -Gu] and 
% form Rf1 = Q1*[Gf;0] and Rw1 = Q1*[Gw;0]
Q1 = [eye(p) -Gu]; Rf1 = Gf; Rw1 = Gw;

% Step 2): determine Q{i} and corresponding Rf{i} and Rw{i}

% initialization
nb = size(S,1);      % number of necessary filters 
Q = cell(nb,1); Rf = cell(nb,1); Rw = cell(nb,1);

opts = struct('sdeg',-3,'smarg',-3,'mindeg',true); 
for i = 1:nb
   % perform Procedure AFD or EFD to compute Q{i}
   indd = (S(i,:) == 0); 
   Qi1 = glnull(ss(Rf1(:,indd)));  
   % initialize Q{i}, Rf{i} and Rw{i}
   Qi = Qi1*Q1; Rfi = Qi1*Rf1; Rwi = Qi1*Rw1; 
   if norm(evalfr(Rwi,rand)) > 0.0001
      % compute the quasi-co-outer-co-inner factorization
      [Rwi,Rwo] = goifac(Rwi,struct('tol',1.e-7)); 
      % update 
      Qi = Rwo\Qi; Rfi = Rwo\Rfi; 
   end   
   % update the solution if [Q{i}, Rf{i}, Rw{i}] is improper or unstable
   [Qi_Rfi_Rwi,M] = glcf([Qi Rfi Rwi],opts);
   % adjust denominator M to unit infinity norm to match example
   Mnorm = norm(M,inf); 
   Q{i} = tf(Qi_Rfi_Rwi(:,1:p+mu)/Mnorm);
   Rf{i} = tf(Qi_Rfi_Rwi(:,p+mu+1:p+mu+mf)/Mnorm);
   Rw{i} = tf(Qi_Rfi_Rwi(:,p+mu+mf+1:end)/Mnorm);
end
Q{1:end}, Rf{1:end}, Rw{1:end}
