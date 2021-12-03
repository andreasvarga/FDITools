% Example 2.2 - Convert a LPV model to LTI form
% Uses the Robust Control Toolbox 

% define rho1 and rho2 as uncertain parameters
r1=ureal('rho1',0); r2=ureal('rho2',0);

% define E, A(rho1,rho2), Bu, C, Du
n = 3; mu = 2; p = 2;       % enter dimensions
E = eye(n);
A = [ -.8 0 0;
      0 -0.5*(1+r1) 0.6*(1+r2);
      0 -0.6*(1+r2) -0.5*(1+r1) ];
Bu = [ 1 1; 1 0; 0 1]; C = [0 1 1;1 1 0]; Du = zeros(p,mu);

% build S(rho)
S = [ E A Bu; zeros(p,n) C Du];

% compute the elements of LFT-based representation
[M,Delta] = lftdata(S);
nd = size(Delta,1);              % size of Delta

% computes orthogonal basis for the range of M21
U1 = orth(M(nd+1:end,1:nd));     %  U1 directly from SVD

% compute Bw and Dw, and define number of noise inputs
Bw = U1(1:n,:); Dw = U1(n+1:end,:); mw = size(U1,2);
